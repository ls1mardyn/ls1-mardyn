/*
 * InSituMegamol.h
 *
 *  Created on: 19 Jul 2018
 *      Author: Oliver Fernandes
 */

#include "InSituMegamol.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

using Log::global_log;

#ifdef INSITU
InSituMegamol::InSituMegamol(void) {
	constexpr bool uint64_tIsSize_t = std::is_same<uint64_t, size_t>::value;
	mardyn_assert(uint64_tIsSize_t);
}

void InSituMegamol::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
	_zmqManager.setConnection(_connectionName);	
	_zmqManager.setReplyBufferSize(_replyBufferSize);
	_zmqManager.setSyncTimeout(_syncTimeout);
	_zmqManager.setModuleNames(domainDecomp->getRank());
	_zmqManager.triggerChainModuleSetup();
}

void InSituMegamol::readXML(XMLfileUnits& xmlconfig) {
	_snapshotInterval = 20;
	xmlconfig.getNodeValue("snapshotInterval", _snapshotInterval);
	global_log->info() << "    ISM: Snapshot interval: "
	        << _snapshotInterval << std::endl;

	_connectionName.assign("tcp://localhost:33333");
	xmlconfig.getNodeValue("connectionName", _connectionName);
	global_log->info() << "    ISM: Connecting to Megamol on: <" 
	        << _connectionName << ">" << std::endl;

	_replyBufferSize = 16384;
	xmlconfig.getNodeValue("replyBufferSize", _replyBufferSize);
	global_log->info() << "    ISM: Megamol reply buffer size (defaults to 16384 byte): "
	        << _replyBufferSize << std::endl;
	_syncTimeout = 10;
	xmlconfig.getNodeValue("syncTimeout", _syncTimeout);
	global_log->info() << "    ISM: Synchronization timeout (s): "
	        << _syncTimeout << std::endl;
}

///reset the particle container with a saved snapshot
void InSituMegamol::beforeEventNewTimestep(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep) {
}

void InSituMegamol::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	if (!(simstep % _snapshotInterval)) {
		//get bbox 
		float bbox[6] {
			0.0f, 0.0f, 0.0f,
			static_cast<float>(domain->getGlobalLength(0)),
			static_cast<float>(domain->getGlobalLength(1)),
			static_cast<float>(domain->getGlobalLength(2))
		};
		//convert simulation time
		float simTime = static_cast<float>(global_simulation->getSimulationTime());
		//build data lists
		std::vector< std::vector<char> > particleLists;
		particleLists.push_back(_buildMmpldDataList(particleContainer));
		//generate seekTable
		std::vector<char> seekTable = _generateMmpldSeekTable(particleLists);

		//dump everything into the _mmpldBuffer
		_resetMmpldBuffer();
		_addMmpldHeader(bbox, simTime);
		_addMmpldSeekTable(seekTable);
		mardyn_assert(*reinterpret_cast<size_t*>(seekTable.data()) == _mmpldBuffer.size());
		_addMmpldFrame(particleLists);
		mardyn_assert(*reinterpret_cast<size_t*>(seekTable.data()+sizeof(size_t)) == _mmpldBuffer.size());
		std::string fname = _writeMmpldBuffer(domainDecomp->getRank());
		_zmqManager.triggerUpdate(fname);
	}
}

std::vector<char> InSituMegamol::_generateMmpldSeekTable(std::vector< std::vector<char> >& particleLists) {
	size_t const headerSize = 60;
	// seek table contains exactly 2 uint64_t, the start of the frame and the end
	size_t frameStart = 60+2*sizeof(size_t);
	// each frame start with an uint32_t, containing the number of particle lists
	// hence, add 4 bytes for that
	size_t frameEnd = frameStart+sizeof(unsigned int);
	// sum the individual particle list counts up
	for (auto const& list : particleLists) {
		frameEnd += list.size();
	}
	std::vector<char> seekTable;
	std::copy(reinterpret_cast<char*>(&frameStart),
			  reinterpret_cast<char*>(&frameStart)+sizeof(size_t),
			  std::back_inserter(seekTable));
	std::copy(reinterpret_cast<char*>(&frameEnd),
			  reinterpret_cast<char*>(&frameEnd)+sizeof(size_t),
			  std::back_inserter(seekTable));
	mardyn_assert(2*sizeof(size_t) == seekTable.size());
	return seekTable;
}

void InSituMegamol::_resetMmpldBuffer(void) {
	_mmpldBuffer.clear();
}

void InSituMegamol::_addMmpldHeader(float bbox[6], float simTime) {
	/// add the standard header data here
	std::string magicId("MMPLD");
	std::copy(magicId.begin(), magicId.end(), std::back_inserter(_mmpldBuffer));
	_mmpldBuffer.push_back(0);
	mardyn_assert(_mmpldBuffer.size() == 6);
	mardyn_assert(sizeof(unsigned short) == 2);
	unsigned short version = 1 * 100 + 3;
	std::copy(reinterpret_cast<char*>(&version),
	          reinterpret_cast<char*>(&version)+sizeof(unsigned short),
			  std::back_inserter(_mmpldBuffer));
	mardyn_assert(_mmpldBuffer.size() == 8);
	mardyn_assert(sizeof(unsigned int) == 4);
	unsigned int numberOfTimesteps = 1;
	std::copy(reinterpret_cast<char*>(&numberOfTimesteps),
	          reinterpret_cast<char*>(&numberOfTimesteps)+sizeof(unsigned int),
			  std::back_inserter(_mmpldBuffer));
	mardyn_assert(_mmpldBuffer.size() == 12);
	std::copy(reinterpret_cast<char*>(bbox),
	          reinterpret_cast<char*>(bbox)+6*sizeof(float),
			  std::back_inserter(_mmpldBuffer));
	mardyn_assert(_mmpldBuffer.size() == 36);
	//add another copy of bbox as fake clipping box
	std::copy(reinterpret_cast<char*>(bbox),
	          reinterpret_cast<char*>(bbox)+6*sizeof(float),
			  std::back_inserter(_mmpldBuffer));
	mardyn_assert(_mmpldBuffer.size() == 60);
}

void InSituMegamol::_addMmpldSeekTable(std::vector<char> seekTable) {
	std::copy(seekTable.begin(), seekTable.end(), std::back_inserter(_mmpldBuffer));
}

void InSituMegamol::_addMmpldFrame(std::vector< std::vector<char> > particleLists) {
	unsigned int numLists = static_cast<unsigned int>(particleLists.size());
	std::copy(reinterpret_cast<char*>(&numLists),
	          reinterpret_cast<char*>(&numLists)+sizeof(float),
			  std::back_inserter(_mmpldBuffer));
	for (auto const& list : particleLists) {
		std::copy(list.begin(), list.end(), std::back_inserter(_mmpldBuffer));
	}
}

std::vector<char> InSituMegamol::_buildMmpldDataList(ParticleContainer* particleContainer) { 
	// add list header
	std::vector<char> dataList;
	dataList.push_back(1); //vertex type
	dataList.push_back(0); //color type
#warning Most particle list meta data is hardcoded (e.g. radius, global color)
	// insert global radius
	float radius = 1.0f;
	std::copy(reinterpret_cast<char*>(&radius),
	          reinterpret_cast<char*>(&radius)+sizeof(float),
			  std::back_inserter(dataList));
	// dummies for Global RGB Color
	dataList.insert(dataList.end(), 4, static_cast<unsigned char>(0.8*255.0));

	// number of particles
	size_t particleCount = static_cast<size_t>(particleContainer->getNumberOfParticles());
	mardyn_assert(sizeof(size_t) == 8);
	std::copy(reinterpret_cast<char*>(&particleCount),
	          reinterpret_cast<char*>(&particleCount)+sizeof(size_t),
			  std::back_inserter(dataList));
	// list bounding box OMG, this needs fixing i guess
	dataList.insert(dataList.end(), 6*sizeof(float), 0);

	// add vertex data
	for (ParticleIterator mol = particleContainer->iterator(); mol.hasNext(); mol.next()) {
		float pos[3] {
			static_cast<float>(mol->r(0)),
			static_cast<float>(mol->r(1)),
			static_cast<float>(mol->r(2))
		};
		std::copy(reinterpret_cast<char*>(&pos),
				  reinterpret_cast<char*>(&pos)+3*sizeof(float),
				  std::back_inserter(dataList));
	}
	return dataList;
}

std::string InSituMegamol::_writeMmpldBuffer(int rank) {
	// all data should be stored, fill in post-data-collection values
	ofstream mmpldFile;
	std::stringstream fname;
	fname << "/dev/shm/part_rnk" << std::setfill('0') << std::setw(6) << rank << ".mmpld";
	mmpldFile.open(fname.str(), std::ios::binary | std::ios::trunc);
	mmpldFile.write(_mmpldBuffer.data(), _mmpldBuffer.size());
	mmpldFile.close();
	global_log->info() << "    ISM: Shared memory file written." << std::endl;
	return fname.str();
}

InSituMegamol::ZmqManager::ZmqManager(void) 
		: _context(zmq_ctx_new(), &zmq_ctx_destroy)
		, _requester(zmq_socket(_context.get(), ZMQ_REQ), &zmq_close) 
		, _publisher(zmq_socket(_context.get(), ZMQ_PUB), &zmq_close) {
	mardyn_assert(_context.get());
	mardyn_assert(_requester.get());
	mardyn_assert(_publisher.get());
	int lingerTime=0;
	zmq_setsockopt(_requester.get(), ZMQ_LINGER, &lingerTime, sizeof(int));
	global_log->info() << "    ISM: Acquired ZmqManager resources." << std::endl;
}

InSituMegamol::ZmqManager::~ZmqManager(void) {
	_publisher.reset();
	_requester.reset();
	_context.reset();
}

InSituMegamol::ZmqManager::ZmqManager(ZmqManager const& rhs) 
		: _context(zmq_ctx_new(), &zmq_ctx_destroy)
		, _requester(zmq_socket(_context.get(), ZMQ_REQ), &zmq_close) 
		, _publisher(zmq_socket(_context.get(), ZMQ_PUB), &zmq_close) {
	
}

void InSituMegamol::ZmqManager::getZmqVersion(void) const {
	global_log->info() << "    ISM: Using ZeroMQ version: " 
	                   << ZMQ_VERSION_MAJOR << "."
	                   << ZMQ_VERSION_MINOR << "."
	                   << ZMQ_VERSION_PATCH << std::endl;
}

void InSituMegamol::ZmqManager::setConnection(std::string connectionName) {
	if (zmq_connect(_requester.get(), connectionName.data())) {
		global_log->info() << "    ISM: Connection failed, releasing resources." << std::endl;
		//TODO stop plugin
	}
	global_log->info() << "    ISM: Connecting to MegaMol on port: " << connectionName << std::endl;
}

void InSituMegamol::ZmqManager::setModuleNames(int rank) {
	_datTag << "::dat" << std::setfill('0') << std::setw(6) << rank;
	_geoTag << "::geo" << std::setfill('0') << std::setw(6) << rank;
}

void InSituMegamol::ZmqManager::triggerChainModuleSetup(void) {
	InSituMegamol::ISM_SYNC_STATUS ismStatus = ISM_SYNC_SYNCHRONIZING;
	do {
		// try to sync to Megamol. Abort either on error or after _syncTimeout sec.
		ismStatus = synchronizeMegamol();
		global_log->info() << "ismStatus: " << ismStatus << std::endl;
		if (ismStatus == ISM_SYNC_TIMEOUT
				|| ismStatus == ISM_SYNC_REPLY_BUFFER_OVERFLOW
				|| ismStatus == ISM_SYNC_UNKNOWN_ERROR) {
			global_log->info() << "    ISM: Megamol sync failed." << std::endl;
			return;
		}
	} while (ismStatus != ISM_SYNC_SUCCESS);

	std::string reply;
	int replySize = 0;
	std::stringstream msg;
	mardyn_assert(_replyBuffer.size() == static_cast<size_t>(_replyBufferSize));
	msg << "mmCreateModule(\"MMPLDDataSource\", \"" << _datTag.str() << "\")\n"
		<< "mmCreateModule(\"OSPRayNHSphereGeometry\", \"" << _geoTag.str() << "\")\n"
		<< "mmCreateCall(\"MultiParticleDataCall\", \"" << _geoTag.str() << "::getData\", \"" << _datTag.str() << "::getData\")\n"
		<< "mmCreateCall(\"CallOSPRayMaterial\", \""<< _geoTag.str() <<"::getMaterialSlot\", " << "\"::mat::deployMaterialSlot\")\n"
		<< "mmCreateChainCall(\"CallOSPRayStructure\", \"::rnd::getStructure\", \"" << _geoTag.str() << "::deployStructureSlot\")\n";
	global_log->info() << "    ISM: Sending creation message\n" << msg.str() << std::endl;
	zmq_send(_requester.get(), msg.str().data(), msg.str().size(), 0);
	replySize = zmq_recv(_requester.get(), _replyBuffer.data(), _replyBufferSize, 0);
	// if (replySize > _replyBufferSize) {
	// 	//some errror
	// }
	// else {
	// 	reply.assign(_replyBuffer.data(), 0, replySize);
	// 	if (replySize == -1 || replySize == 0) {
	// 		global_log->info() << "    ISM: ZMQ reply was empty.("<< replySize <<")" << std::endl;
	// 	}
	// 	else {
	// 		global_log->info() << "    ISM: ZMQ reply: " << replySize << std::endl;
	// 	}
	// }
}

InSituMegamol::ISM_SYNC_STATUS InSituMegamol::ZmqManager::synchronizeMegamol(void) {
	std::string msg("return mmListModules()");
	std::string reply;
	global_log->info() << "    ISM: Requesting module list" << std::endl;
	zmq_send(_requester.get(), msg.data(), msg.size(), ZMQ_DONTWAIT);
	int replySize = zmq_recv(_requester.get(), _replyBuffer.data(), _replyBufferSize, ZMQ_DONTWAIT);
	InSituMegamol::ISM_SYNC_STATUS currentStatus;
	// evaluate Megamol's reply
	if (replySize > _replyBufferSize) {
		// This is an error. Return and handle outside. Should probably do an exception here.
		global_log->info() << "    ISM: reply size exceeded buffer size." << std::endl;
		return InSituMegamol::ISM_SYNC_REPLY_BUFFER_OVERFLOW;
	}
	if (replySize == -1) {
		// No reply from Megamol yet, let's wait.
		global_log->info() << "    ISM: Megamol not ready.("<< replySize <<")" << std::endl;
		currentStatus = InSituMegamol::ISM_SYNC_SYNCHRONIZING;
		if ((_timeoutClock()) > _syncTimeout) {
			// we waited long enough, consider the sync failed
			return InSituMegamol::ISM_SYNC_TIMEOUT;
		}
		else {
			return InSituMegamol::ISM_SYNC_SYNCHRONIZING;
		}
	} 
	if (replySize == 0) {
		// Megamol replied, but does not seem ready yet, let's wait.
		global_log->info() << "    ISM: ZMQ reply was empty.("<< replySize <<")" << std::endl;
		if (_timeoutClock() > _syncTimeout) {
			// we waited long enough, consider the sync failed
			return InSituMegamol::ISM_SYNC_TIMEOUT;
		}
		else {
			return InSituMegamol::ISM_SYNC_SYNCHRONIZING;
		}
	}
	if (replySize > 0 && replySize < _replyBufferSize) {
		// Check Megamol's reply
		reply.assign(_replyBuffer.data(), 0, replySize);
		global_log->info() << "    ISM: ZMQ reply received." << std::endl;
		if (reply.find("OSPRayRenderer") != std::string::npos) {
			// Megamol seems ready. Confirm synchronization.
			global_log->info() << "    ISM: Synchronized Megamol." << std::endl;
			return InSituMegamol::ISM_SYNC_SUCCESS;
		}
		else {
 			// Megamol replied, but does not seem ready yet, let's wait.
			if (_timeoutClock() > _syncTimeout) {
				// we waited long enough, consider the sync failed
				return InSituMegamol::ISM_SYNC_TIMEOUT;
			}
			else {
				return InSituMegamol::ISM_SYNC_SYNCHRONIZING;
			}
		}
	}
	return InSituMegamol::ISM_SYNC_UNKNOWN_ERROR;
}

void InSituMegamol::ZmqManager::triggerUpdate(std::string fname) {
	std::stringstream msg;
	msg << "mmSetParamValue(\"" << _datTag.str() << "::filename\", \"" << fname << "\")";
	zmq_send(_requester.get(), msg.str().data(), msg.str().size(), ZMQ_DONTWAIT);
	zmq_recv(_requester.get(), _replyBuffer.data(), _replyBufferSize, ZMQ_DONTWAIT);
}

#else

InSituMegamol::InSituMegamol(void) {
	global_log->info() << "InSituMegamol: This is a just a dummy."
			<< "Set INSITU=1 on make command line to enable the actual plugin."
			<< std::endl;
}

#endif //INSITU