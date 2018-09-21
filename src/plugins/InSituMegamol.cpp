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
#define BLOCK_POLICY_HANDSHAKE 0
#define BLOCK_POLICY_UPDATE ZMQ_DONTWAIT

InSituMegamol::InSituMegamol(void) 
		: _isEnabled(false) {
	constexpr bool uint64_tIsSize_t = std::is_same<uint64_t, size_t>::value;
	mardyn_assert(uint64_tIsSize_t);
}

void InSituMegamol::init(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
	_zmqManager.setConnection(_connectionName);	
	_zmqManager.setReplyBufferSize(_replyBufferSize);
	_zmqManager.setSyncTimeout(_syncTimeout);
	_zmqManager.setModuleNames(domainDecomp->getRank());
	_createFnameRingBuffer(domainDecomp->getRank(), _ringBufferSize);
	_isEnabled = _zmqManager.performHandshake();
	std::stringstream llFname;
	llFname << "/dev/shm/local" << std::setfill('0') << std::setw(6) << ".log";
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
	_ringBufferSize = 5;
	xmlconfig.getNodeValue("ringBufferSize", _ringBufferSize);
	global_log->info() << "    ISM: Ring buffer size: "
	        << _ringBufferSize << std::endl;
}

///reset the particle container with a saved snapshot
void InSituMegamol::beforeEventNewTimestep(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep) {
}

void InSituMegamol::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	if (!(simstep % _snapshotInterval)) {
		if (!_isEnabled) {
			global_log->info() << "    ISM: Disabled. Skipping InSitu plugin." << std::endl;
			return;
		}
		auto start = std::chrono::high_resolution_clock::now();
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
		std::string fname = _writeMmpldBuffer(domainDecomp->getRank(), simstep);
		auto end = std::chrono::high_resolution_clock::now();
		_addTimerEntry(simstep, std::chrono::duration<double>(end-start).count());
		_isEnabled = _zmqManager.triggerUpdate(fname);
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
	// dummies for Global RGB Color, insert 0.8*255 4 times for RGBA.
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
	for (auto mol = particleContainer->iterator(); mol.isValid(); ++mol) {
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

std::string InSituMegamol::_writeMmpldBuffer(int rank, unsigned long simstep) {
	// all data should be stored, fill in post-data-collection values
	ofstream mmpldFile;
	// std::stringstream fname;
	// fname << "/dev/shm/part_rnkI" << std::setfill('0') << std::setw(6) << rank 
	// 		<< "_T" << std::setfill('0') << std::setw(7) << simstep << ".mmpld";
	std::string fname(_getNextFname());
	mmpldFile.open(fname, std::ios::binary | std::ios::trunc);
	mmpldFile.write(_mmpldBuffer.data(), _mmpldBuffer.size());
	mmpldFile.close();
	global_log->info() << "    ISM: Shared memory file written." << std::endl;
	return fname;
}

void InSituMegamol::_createFnameRingBuffer(int const rank, int const size) {
	std::stringstream fname;
	for (size_t i=0; i<5; ++i) {
		fname << "/dev/shm/part_rnk" << std::setfill('0') << std::setw(6) << rank 
				<< "_buf" << std::setfill('0') << std::setw(2) << i << ".mmpld";
		_fnameRingBuffer.push_back(fname.str());
		fname.str("");
	}
}

std::string InSituMegamol::_getNextFname(void) {
	static RingBuffer::iterator nextFname = _fnameRingBuffer.begin();
	auto temp = nextFname;
	if (++nextFname == _fnameRingBuffer.end()) {
		nextFname = _fnameRingBuffer.begin();
	}
	return *temp;
}

void InSituMegamol::_addTimerEntry(unsigned long simstep, double secs) {
	std::ofstream localLog;
	std::stringstream timerLine;
	// dump times in microseconds
	timerLine << "T: " << std::setfill(' ') << std::setw(8) <<  simstep 
			<< " d: " << std::fixed << std::setprecision(6);
	localLog.open(_localLogFname, std::ios::ate);
	localLog << timerLine.str() << "\n";
	localLog.close();
}

InSituMegamol::ZmqManager::ZmqManager(void) 
		: _sendCount(0)
		, _context(zmq_ctx_new(), &zmq_ctx_destroy)
		, _requester(zmq_socket(_context.get(), ZMQ_REQ), &zmq_close) 
		, _publisher(zmq_socket(_context.get(), ZMQ_PUB), &zmq_close) {
	mardyn_assert(_context.get());
	mardyn_assert(_requester.get());
	mardyn_assert(_publisher.get());
	int lingerTime=0;
	zmq_setsockopt(_requester.get(), ZMQ_LINGER, &lingerTime, sizeof(int));
	global_log->info() << "    ISM: Acquired ZmqManager resources." << std::endl;
	global_log->info() << "    ISM: Using ZeroMQ version: " << getZmqVersion() << std::endl;
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

std::string InSituMegamol::ZmqManager::getZmqVersion(void) const {
	std::stringstream version;
	version << ZMQ_VERSION_MAJOR << "."
	        << ZMQ_VERSION_MINOR << "."
	        << ZMQ_VERSION_PATCH;
	return version.str();
}

void InSituMegamol::ZmqManager::setConnection(std::string connectionName) {
	if (zmq_connect(_requester.get(), connectionName.data())) {
		global_log->info() << "    ISM: Connection failed, releasing resources." << std::endl;
		// TODO: propagate plugin shutdown;
	}
	global_log->info() << "    ISM: Requester initialized." << std::endl;
}

void InSituMegamol::ZmqManager::setModuleNames(int rank) {
	_datTag << "::dat" << std::setfill('0') << std::setw(6) << rank;
	_geoTag << "::geo" << std::setfill('0') << std::setw(6) << rank;
}

bool InSituMegamol::ZmqManager::performHandshake(void) {
	InSituMegamol::ISM_SYNC_STATUS ismStatus = ISM_SYNC_SYNCHRONIZING;
	do {
		// try to sync to Megamol. Abort either on error or after _syncTimeout sec.
		ismStatus = _synchronizeMegamol();
		if (ismStatus == ISM_SYNC_TIMEOUT
				|| ismStatus == ISM_SYNC_REPLY_BUFFER_OVERFLOW
				|| ismStatus == ISM_SYNC_SEND_FAILED
				|| ismStatus == ISM_SYNC_UNKNOWN_ERROR) {
			global_log->info() << "    ISM: Megamol sync failed." << std::endl;
			return false;
		}
	} while (ismStatus != ISM_SYNC_SUCCESS);

	std::string reply;
	std::stringstream msg;
	mardyn_assert(_replyBuffer.size() == static_cast<size_t>(_replyBufferSize));
	msg << "mmCreateModule(\"MMPLDDataSource\", \"" << _datTag.str() << "\")\n"
		<< "mmCreateModule(\"OSPRayNHSphereGeometry\", \"" << _geoTag.str() << "\")\n"
		<< "mmCreateCall(\"MultiParticleDataCall\", \"" << _geoTag.str() << "::getData\", \"" << _datTag.str() << "::getData\")\n"
		<< "mmCreateCall(\"CallOSPRayMaterial\", \""<< _geoTag.str() <<"::getMaterialSlot\", " << "\"::mat::deployMaterialSlot\")\n"
		<< "mmCreateChainCall(\"CallOSPRayStructure\", \"::rnd::getStructure\", \"" << _geoTag.str() << "::deployStructureSlot\")\n";
	global_log->info() << "    ISM: Sending creation message." << std::endl;
	zmq_send(_requester.get(), msg.str().data(), msg.str().size(), BLOCK_POLICY_HANDSHAKE);
	int replySize = zmq_recv(_requester.get(), _replyBuffer.data(), _replyBufferSize, BLOCK_POLICY_HANDSHAKE);
	// surprisingly, stuff actually worked, enable the plugin
	return true;
}

InSituMegamol::ISM_SYNC_STATUS InSituMegamol::ZmqManager::_synchronizeMegamol(void) {
	std::string msg("return mmListModules()");
	std::string reply;
	std::stringstream statusStr("");
	statusStr << "    ISM: Requesting module list...";
	if (zmq_send(_requester.get(), msg.data(), msg.size(), BLOCK_POLICY_HANDSHAKE) == -1) {
		statusStr << "send message failed. Error: " << strerror(errno);
		global_log->info() << statusStr.str() << std::endl;
		return InSituMegamol::ISM_SYNC_SEND_FAILED;
	}
	int replySize = zmq_recv(_requester.get(), _replyBuffer.data(), _replyBufferSize, BLOCK_POLICY_HANDSHAKE);
	// evaluate Megamol's reply
	if (replySize == -1) {
		// No reply from Megamol yet, let's wait.
		statusStr << "Megamol not responding. Error: " << strerror(errno);
		global_log->info() << statusStr.str() << std::endl;
		return _timeoutCheck();
	} 
	if (replySize > _replyBufferSize) {
		// This is an error. Return and handle outside. Should probably do an exception here.
		global_log->info() << "    ISM: reply size exceeded buffer size." << std::endl;
		return InSituMegamol::ISM_SYNC_REPLY_BUFFER_OVERFLOW;
	}
	if (replySize == 0) {
		// Megamol replied, but does not seem ready yet, let's wait.
		statusStr << "ZMQ reply was empty.";
		global_log->info() << statusStr.str() << std::endl;
		return _timeoutCheck();
	}
	if (replySize > 0 && replySize < _replyBufferSize) {
		// Check Megamol's reply
		reply.assign(_replyBuffer.data(), 0, replySize);
		statusStr << "ZMQ reply received (size: " << replySize << ").";
		global_log->info() << statusStr.str() << std::endl;
		if (reply.find("OSPRayRenderer") != std::string::npos) {
			// Megamol seems ready. Confirm synchronization happened.
			global_log->info() << "    ISM: Synchronized Megamol." << std::endl;
			return InSituMegamol::ISM_SYNC_SUCCESS;
		}
		else {
 			// Megamol replied, but does not seem ready yet, let's wait.
			return _timeoutCheck();
		}
	}
	return InSituMegamol::ISM_SYNC_UNKNOWN_ERROR;
}

InSituMegamol::ISM_SYNC_STATUS InSituMegamol::ZmqManager::_timeoutCheck(void) const {
	static int iterationCounter = 0;
	std::chrono::high_resolution_clock hrc;
	auto again = hrc.now() + std::chrono::milliseconds(1000);
	while (hrc.now() < again);
	++iterationCounter;
	// return iterationCounter > _syncTimeout;
	if (iterationCounter > _syncTimeout) {
		// we waited long enough, consider the sync failed
		return InSituMegamol::ISM_SYNC_TIMEOUT;
	}
	else {
		return InSituMegamol::ISM_SYNC_SYNCHRONIZING;
	}
}

bool InSituMegamol::ZmqManager::triggerUpdate(std::string fname) {
	std::stringstream msg;
	msg << "mmSetParamValue(\"" << _datTag.str() << "::filename\", \"" << fname << "\")";
	_send(msg.str(), BLOCK_POLICY_HANDSHAKE);
	int replySize = _recv(BLOCK_POLICY_UPDATE);
	std::stringstream statusStr;
	if (replySize > _replyBufferSize) {
		// This is an error. Should probably throw here.
		global_log->info() << "    ISM: reply size exceeded buffer size. Disabling plugin." << std::endl;
		return false;
	}
	if (replySize == -1) {
		// Reply failed, raise unhandled reply counter
		statusStr << "    ISM: Megamol not ready (size: -1). ";
		global_log->info() << statusStr.str() << std::endl;
		return true;
	} 
	if (replySize == 0) {
		// Megamol replied with 0. This is what we want.
		statusStr << "    ISM: ZMQ reply was empty (size: 0).";
		global_log->info() << statusStr.str() << std::endl;
		return true;
	}
	if (replySize > 0 && replySize < _replyBufferSize) {
		// Megamol's reply was not empty. Dump the reply, stop the plugin.
		std::string reply;
		reply.assign(_replyBuffer.data(), 0, replySize);
		statusStr << "    ISM: ZMQ reply received (size: " << replySize << ", should be 0).";
		statusStr << "    ISM: Reply: <" << reply << ">. Disabling plugin.";
		global_log->info() << statusStr.str() << std::endl;
		return false;
	}
	return false;
}

int InSituMegamol::ZmqManager::_send(std::string msg, int blockPolicy) {
	int status;
	global_log->info() << "    ISM: _sendCount: " << _sendCount << std::endl;
	while (_sendCount > 0) {
		status = _recv(blockPolicy);
		if (status == -1 && errno != EAGAIN) {
			global_log->info() << "    ISM: Stuff is really messed up." << std::endl;
			return status;
		}
	}
	status = zmq_send(_requester.get(), msg.data(), msg.size(), blockPolicy);
	if (status != -1) {
		++_sendCount;
	}
	return status;
}

int InSituMegamol::ZmqManager::_recv(int blockPolicy) {
	int status = zmq_recv(_requester.get(), _replyBuffer.data(), _replyBufferSize, blockPolicy);
	if (status != -1) {
		--_sendCount;
	}
	return status;
}
#else

InSituMegamol::InSituMegamol(void) {
	global_log->info() << "InSituMegamol: This is a just a dummy."
			<< "Set INSITU=1 on make command line to enable the actual plugin."
			<< std::endl;
}

#endif //INSITU