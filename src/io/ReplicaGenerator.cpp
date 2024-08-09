#include "ReplicaGenerator.h"

/* the following macro has to be defined to use math constants in cmath */
#define _USE_MATH_DEFINES  1

#include <cmath>
#include <vector>
#include <cstdint>
#include <numeric>

#include "Domain.h"
#include "Simulation.h"

#include "parallel/DomainDecompBase.h"

#ifdef ENABLE_MPI
#include "parallel/ParticleData.h"
#include "parallel/DomainDecomposition.h"
#endif

#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/mardyn_assert.h"
#include "utils/Logger.h"


enum MoleculeFormat : uint32_t {
	ICRVQD, IRV, ICRV
};

ReplicaGenerator::ReplicaGenerator()
		:
		_numParticlesTotal(0),
		_numBlocksXZ(0),
		_numBlocksLiqY(0),
		_numBlocksVapY(0),
		_nIndexLiqBeginY(0),
		_nIndexLiqEndY(0),
		_nMoleculeFormat(ICRVQD),
		_moleculeDataReader(nullptr),
		_dMoleculeDiameter(0.),
		_fspY{0., 0., 0., 0., 0., 0.},
		_nSystemType(ST_UNKNOWN) {
	// init component ID change data structures
	uint32_t numComponents = global_simulation->getEnsemble()->getComponents()->size();
	_vecChangeCompIDsVap.resize(numComponents);
	_vecChangeCompIDsLiq.resize(numComponents);
	std::iota(std::begin(_vecChangeCompIDsVap), std::end(_vecChangeCompIDsVap), 0);
	std::iota(std::begin(_vecChangeCompIDsLiq), std::end(_vecChangeCompIDsLiq), 0);
}

ReplicaGenerator::~ReplicaGenerator() = default;

void ReplicaGenerator::readReplicaPhaseSpaceHeader(SubDomain& subDomain) {
	XMLfileUnits inp(subDomain.strFilePathHeader);

	if(not inp.changecurrentnode("/mardyn")) {
		Log::global_log->error() << "Could not find root node /mardyn in XML input file." << std::endl;
		Log::global_log->fatal() << "Not a valid MarDyn XML input file." << std::endl;
		mardyn_exit(1);
	}

	bool bInputOk = true;
	double dCurrentTime = 0.;
	double dBL[3];
	std::string strMoleculeFormat;
	bInputOk = bInputOk && inp.changecurrentnode("headerinfo");
	bInputOk = bInputOk && inp.getNodeValue("time", dCurrentTime);
	bInputOk = bInputOk && inp.getNodeValue("length/x", dBL[0]);
	bInputOk = bInputOk && inp.getNodeValue("length/y", dBL[1]);
	bInputOk = bInputOk && inp.getNodeValue("length/z", dBL[2]);
	bInputOk = bInputOk && inp.getNodeValue("number", subDomain.numParticles);
	bInputOk = bInputOk && inp.getNodeValue("format@type", strMoleculeFormat);
	subDomain.dVolume = 1;
	for(uint8_t di = 0; di < 3; ++di) {
		subDomain.arrBoxLength.at(di) = dBL[di];
		subDomain.dVolume *= dBL[di];
	}
	subDomain.dDensity = subDomain.numParticles / subDomain.dVolume;

	if(not bInputOk) {
		Log::global_log->error() << "Content of file: '" << subDomain.strFilePathHeader << "' corrupted! Program exit ..."
							<< std::endl;
		mardyn_exit(1);
	}

	if("ICRVQD" == strMoleculeFormat)
		_nMoleculeFormat = ICRVQD;
	else if("IRV" == strMoleculeFormat)
		_nMoleculeFormat = IRV;
	else if("ICRV" == strMoleculeFormat)
		_nMoleculeFormat = ICRV;
	else {
		Log::global_log->error() << "Not a valid molecule format: " << strMoleculeFormat << ", program exit ..." << std::endl;
		mardyn_exit(1);
	}
}

void ReplicaGenerator::readReplicaPhaseSpaceData(SubDomain& subDomain, DomainDecompBase* domainDecomp) {
#ifdef ENABLE_MPI
	if(domainDecomp->getRank() == 0) {
#endif
	subDomain.strFilePathData = string_utils::trim(subDomain.strFilePathData);
	Log::global_log->info() << "Opening phase space file " << subDomain.strFilePathData << std::endl;
	std::ifstream ifs;
	ifs.open(subDomain.strFilePathData.c_str(), std::ios::binary | std::ios::in);
	if(!ifs.is_open()) {
		Log::global_log->error() << "Could not open phaseSpaceFile " << subDomain.strFilePathData << std::endl;
		mardyn_exit(1);
	}

	Log::global_log->info() << "Reading phase space file " << subDomain.strFilePathData << std::endl;

	std::vector<Component>& components = *(_simulation.getEnsemble()->getComponents());

	// Select appropriate reader
	switch (_nMoleculeFormat) {
		case ICRVQD:
			_moleculeDataReader = new MoleculeDataReaderICRVQD();
			break;
		case ICRV:
			_moleculeDataReader = new MoleculeDataReaderICRV();
			break;
		case IRV:
			_moleculeDataReader = new MoleculeDataReaderIRV();
			break;
	}

	for(uint64_t pi = 0; pi < subDomain.numParticles; pi++) {
		Molecule mol;
		_moleculeDataReader->read(ifs, mol, components);
		subDomain.vecParticles.push_back(mol);
	}
#ifdef ENABLE_MPI
	}
#endif

	/* distribute molecules to other MPI processes */
#ifdef ENABLE_MPI
	unsigned long num_particles = subDomain.vecParticles.size();
	MPI_CHECK( MPI_Bcast(&num_particles, 1, MPI_UNSIGNED_LONG, 0, domainDecomp->getCommunicator()) );

#define PARTICLE_BUFFER_SIZE  (16*1024)
	ParticleData particle_buff[PARTICLE_BUFFER_SIZE];
	int particle_buff_pos = 0;
	MPI_Datatype mpi_Particle;
	ParticleData::getMPIType(mpi_Particle);

	if(domainDecomp->getRank() == 0) {
		for(unsigned long i = 0; i < num_particles; ++i) {
			ParticleData::MoleculeToParticleData(particle_buff[particle_buff_pos], subDomain.vecParticles[i]);
			particle_buff_pos++;
			if ((particle_buff_pos >= PARTICLE_BUFFER_SIZE) || (i == num_particles - 1)) {
				Log::global_log->debug() << "broadcasting(sending) particles" << std::endl;
				MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, 0, domainDecomp->getCommunicator());
				particle_buff_pos = 0;
			}
		}
	} else {
		for(unsigned long i = 0; i < num_particles; ++i) {
			if(i % PARTICLE_BUFFER_SIZE == 0) {
				Log::global_log->debug() << "broadcasting(receiving) particles" << std::endl;
				MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, 0, domainDecomp->getCommunicator());
				particle_buff_pos = 0;
			}
			Molecule m;
			ParticleData::ParticleDataToMolecule(particle_buff[particle_buff_pos], m);
			particle_buff_pos++;
			subDomain.vecParticles.push_back(m);
		}
	}
	Log::global_log->debug() << "broadcasting(sending/receiving) particles complete" << std::endl;
#endif
	Log::global_log->info() << "Reading Molecules done" << std::endl;
}

void ReplicaGenerator::readXML(XMLfileUnits& xmlconfig) {
	Log::global_log->debug() << "Reading config for ReplicaGenerator" << std::endl;

	_nSystemType = ST_UNKNOWN;
	std::string strType = "unknown";
	xmlconfig.getNodeValue("type", strType);
	if("homogeneous" == strType) {
		_nSystemType = ST_HOMOGENEOUS;
	} else if("heterogeneous_VLV" == strType) {
		_nSystemType = ST_HETEROGENEOUS_VAPOR_LIQUID_VAPOR;
	} else if("heterogeneous_LV" == strType) {
		_nSystemType = ST_HETEROGENEOUS_LIQUID_VAPOR;
	} else {
		Log::global_log->error() << "Specified wrong type at XML path: " << xmlconfig.getcurrentnodepath() << "/type"
							<< std::endl;
		mardyn_exit(-1);
	}

	SubDomain sd;
	xmlconfig.getNodeValue("files/vapor/header", sd.strFilePathHeader);
	xmlconfig.getNodeValue("files/vapor/data", sd.strFilePathData);
	_vecSubDomains.push_back(sd);
	if(_nSystemType == ST_HETEROGENEOUS_VAPOR_LIQUID_VAPOR || _nSystemType == ST_HETEROGENEOUS_LIQUID_VAPOR) {
		SubDomain sd2;
		xmlconfig.getNodeValue("files/liquid/header", sd2.strFilePathHeader);
		xmlconfig.getNodeValue("files/liquid/data", sd2.strFilePathData);
		_vecSubDomains.push_back(sd2);
	}

	xmlconfig.getNodeValue("numblocks/xz", _numBlocksXZ);
	xmlconfig.getNodeValue("numblocks/vapor", _numBlocksVapY);
	if(_nSystemType == ST_HETEROGENEOUS_VAPOR_LIQUID_VAPOR || _nSystemType == ST_HETEROGENEOUS_LIQUID_VAPOR)
		xmlconfig.getNodeValue("numblocks/liquid", _numBlocksLiqY);

	Log::global_log->info() << "Replicating " << _numBlocksXZ << " x " << _numBlocksXZ << " boxes in XZ layers."
					   << std::endl;

	if(_nSystemType == ST_HETEROGENEOUS_VAPOR_LIQUID_VAPOR) {
		// liquid blocks begin/end index
		_nIndexLiqBeginY = _numBlocksVapY;
		_nIndexLiqEndY = _numBlocksVapY + _numBlocksLiqY - 1;

		xmlconfig.getNodeValue("diameter", _dMoleculeDiameter);
		Log::global_log->info() << "Using molecule diameter: " << _dMoleculeDiameter
						   << " for spacing between liquid and vapour phase. " << std::endl;
	} else if(_nSystemType == ST_HETEROGENEOUS_LIQUID_VAPOR) {
		// liquid blocks begin/end index
		_nIndexLiqBeginY = 0;
		_nIndexLiqEndY = _numBlocksLiqY - 1;

		xmlconfig.getNodeValue("diameter", _dMoleculeDiameter);
		Log::global_log->info() << "Using molecule diameter: " << _dMoleculeDiameter
						   << " for spacing between liquid and vapour phase. " << std::endl;
	}

	// change identity of molecules by component ID (zero based))
	{
		// vapor system
		std::string oldpath = xmlconfig.getcurrentnodepath();
		if(xmlconfig.changecurrentnode("componentIDs/vapor")) {
			uint8_t numChanges = 0;
			XMLfile::Query query = xmlconfig.query("change");
			numChanges = query.card();
			Log::global_log->info() << "Number of components to change: " << (uint32_t) numChanges << std::endl;
			if(numChanges < 1) {
				Log::global_log->error() << "No component change defined in XML-config file. Program exit ..." << std::endl;
				mardyn_exit(-1);
			}
			XMLfile::Query::const_iterator changeIter;
			for(changeIter = query.begin(); changeIter; changeIter++) {
				xmlconfig.changecurrentnode(changeIter);
				uint32_t nFrom, nTo;
				nFrom = nTo = 1;
				xmlconfig.getNodeValue("from", nFrom);
				xmlconfig.getNodeValue("to", nTo);
				_vecChangeCompIDsVap.at(nFrom - 1) = nTo - 1;
			}
		}
		xmlconfig.changecurrentnode(oldpath);

		// liquid system
		oldpath = xmlconfig.getcurrentnodepath();
		if(xmlconfig.changecurrentnode("componentIDs/liquid")) {
			uint8_t numChanges = 0;
			XMLfile::Query query = xmlconfig.query("change");
			numChanges = query.card();
			Log::global_log->info() << "Number of components to change: " << (uint32_t) numChanges << std::endl;
			if(numChanges < 1) {
				Log::global_log->error() << "No component change defined in XML-config file. Program exit ..." << std::endl;
				mardyn_exit(-1);
			}
			XMLfile::Query::const_iterator changeIter;
			for(changeIter = query.begin(); changeIter; changeIter++) {
				xmlconfig.changecurrentnode(changeIter);
				uint32_t nFrom, nTo;
				nFrom = nTo = 1;
				xmlconfig.getNodeValue("from", nFrom);
				xmlconfig.getNodeValue("to", nTo);
				_vecChangeCompIDsLiq.at(nFrom - 1) = nTo - 1;
			}
		}
		xmlconfig.changecurrentnode(oldpath);
	}

	this->init();
}

void ReplicaGenerator::init() {
	DomainDecompBase* domainDecomp = &global_simulation->domainDecomposition();
	Log::global_log->info() << domainDecomp->getRank() << ": Init Replica VLE ..." << std::endl;

	for(auto&& sd : _vecSubDomains) {
		this->readReplicaPhaseSpaceHeader(sd);
		this->readReplicaPhaseSpaceData(sd, domainDecomp);
	}

	// total number of particles
	switch (_nSystemType) {
		case ST_HOMOGENEOUS:
			_numParticlesTotal = _numBlocksVapY * _vecSubDomains.at(0).numParticles * _numBlocksXZ * _numBlocksXZ;
			break;
		case ST_HETEROGENEOUS_VAPOR_LIQUID_VAPOR:
			_numParticlesTotal = (2 * _numBlocksVapY * _vecSubDomains.at(0).numParticles +
								  _numBlocksLiqY * _vecSubDomains.at(1).numParticles) * _numBlocksXZ * _numBlocksXZ;
			break;
		case ST_HETEROGENEOUS_LIQUID_VAPOR:
			_numParticlesTotal = (_numBlocksVapY * _vecSubDomains.at(0).numParticles +
								  _numBlocksLiqY * _vecSubDomains.at(1).numParticles) * _numBlocksXZ * _numBlocksXZ;
			break;
		default:
			_numParticlesTotal = 0;
	}

	// update global length of domain
	double dLength[3];
	switch (_nSystemType) {
		case ST_HOMOGENEOUS:
			dLength[1] = _numBlocksVapY * _vecSubDomains.at(0).arrBoxLength.at(1);
			break;
		case ST_HETEROGENEOUS_VAPOR_LIQUID_VAPOR:
			dLength[1] = 2 * _numBlocksVapY * _vecSubDomains.at(0).arrBoxLength.at(1) +
						 _numBlocksLiqY * _vecSubDomains.at(1).arrBoxLength.at(1);
			break;
		case ST_HETEROGENEOUS_LIQUID_VAPOR:
			dLength[1] = _numBlocksVapY * _vecSubDomains.at(0).arrBoxLength.at(1) +
						 _numBlocksLiqY * _vecSubDomains.at(1).arrBoxLength.at(1);
			break;
		default:
			dLength[1] = 0;
			break;
	}
	dLength[0] = _numBlocksXZ * _vecSubDomains.at(0).arrBoxLength.at(0);
	dLength[2] = _numBlocksXZ * _vecSubDomains.at(0).arrBoxLength.at(2);
	for(uint8_t di = 0; di < 3; ++di)
		global_simulation->getDomain()->setGlobalLength(di, dLength[di]);
	Log::global_log->info() << "Domain box length = " << dLength[0] << ", " << dLength[1] << ", " << dLength[2] << std::endl;

/*
	// Reset domain decomposition
	if (domainDecomp != nullptr) {
		delete domainDecomp;
	}
#ifndef ENABLE_MPI
	Log::global_log->info() << "Initializing the alibi domain decomposition ... " << std::endl;
	domainDecomp = new DomainDecompBase();
#else
	Log::global_log->info() << "Initializing the standard domain decomposition ... " << std::endl;
	domainDecomp = (DomainDecompBase*) new DomainDecomposition();
#endif
	Log::global_log->info() << "Initialization done" << std::endl;
	domainDecomp->readXML(xmlconfig);
*/

	double dPhaseLengthVapY = 0.;
	double dPhaseLengthLiqY = 0.;
	switch (_nSystemType) {
		case ST_HETEROGENEOUS_VAPOR_LIQUID_VAPOR:
			// calc free space positions
			dPhaseLengthVapY = _numBlocksVapY * _vecSubDomains.at(0).arrBoxLength.at(1);
			dPhaseLengthLiqY = _numBlocksLiqY * _vecSubDomains.at(1).arrBoxLength.at(1);
			_fspY[0] = dPhaseLengthVapY - _dMoleculeDiameter;
			_fspY[1] = dPhaseLengthVapY;
			_fspY[2] = dPhaseLengthVapY + dPhaseLengthLiqY;
			_fspY[3] = dPhaseLengthVapY + dPhaseLengthLiqY + _dMoleculeDiameter;
			_fspY[4] = dPhaseLengthVapY + dPhaseLengthLiqY + dPhaseLengthVapY - _dMoleculeDiameter;
			_fspY[5] = dPhaseLengthVapY + dPhaseLengthLiqY + dPhaseLengthVapY;
			break;
		case ST_HETEROGENEOUS_LIQUID_VAPOR:
			// calc free space positions
			dPhaseLengthVapY = _numBlocksVapY * _vecSubDomains.at(0).arrBoxLength.at(1);
			dPhaseLengthLiqY = _numBlocksLiqY * _vecSubDomains.at(1).arrBoxLength.at(1);
			_fspY[0] = 1.;
			_fspY[1] = 0.;
			_fspY[2] = dPhaseLengthLiqY;
			_fspY[3] = dPhaseLengthLiqY + _dMoleculeDiameter;
			_fspY[4] = dPhaseLengthLiqY + dPhaseLengthVapY - _dMoleculeDiameter;
			_fspY[5] = dPhaseLengthLiqY + dPhaseLengthVapY;
			break;
	}
}

unsigned long
ReplicaGenerator::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {
//	if(3 != domainDecomp->getRank() )
//		return 0;
	global_simulation->timers()->start("REPLICA_GENERATOR_VLE_INPUT");
	Log::global_log->info() << "Constructing Replica VLE" << std::endl;

	double bbMin[3];
	double bbMax[3];
	double bbLength[3];
	uint64_t numBlocks[3];
	uint64_t startIndex[3];

	domainDecomp->getBoundingBoxMinMax(domain, bbMin, bbMax);
	for(int d = 0; d < 3; ++d) {
		bbLength[d] = bbMax[d] - bbMin[d];
		numBlocks[d] = ceil(bbLength[d] / _vecSubDomains.at(0).arrBoxLength.at(d) +
							1);  // +1: Particles were missing in some processes without +1 (rounding error??)
		startIndex[d] = floor(bbMin[d] / _vecSubDomains.at(0).arrBoxLength.at(d));
	}

	// Init maxID
	int nRank = domainDecomp->getRank();
	double dVolumeSubdomain = bbLength[0] * bbLength[1] * bbLength[2];
	double dDensityMax = 0;
	for(const auto& sd : _vecSubDomains) {
		if(sd.dDensity > dDensityMax)
			dDensityMax = sd.dDensity;
	}

#ifndef NDEBUG
	std::cout << domainDecomp->getRank() << ": bbMin = " << bbMin[0] << ", " << bbMin[1] << ", " << bbMin[2] << std::endl;
	std::cout << domainDecomp->getRank() << ": bbMax = " << bbMax[0] << ", " << bbMax[1] << ", " << bbMax[2] << std::endl;
	std::cout << domainDecomp->getRank() << ": bbLength = " << bbLength[0] << ", " << bbLength[1] << ", " << bbLength[2]
		 << std::endl;
	std::cout << domainDecomp->getRank() << ": numBlocks = " << numBlocks[0] << ", " << numBlocks[1] << ", " << numBlocks[2]
		 << std::endl;
	std::cout << domainDecomp->getRank() << ": startIndex = " << startIndex[0] << ", " << startIndex[1] << ", "
		 << startIndex[2] << std::endl;
	std::cout << domainDecomp->getRank() << ": BoxLength = " << _vecSubDomains.at(0).arrBoxLength.at(0) << ","
																									  " "
		 << _vecSubDomains.at(0).arrBoxLength.at(1) << ","
													   " " << _vecSubDomains.at(0).arrBoxLength.at(2) << std::endl;
	std::cout << domainDecomp->getRank() << ": bbLength/BoxLength = "
		 << bbLength[0] / _vecSubDomains.at(0).arrBoxLength.at(0) << ","
																	 " "
		 << bbLength[1] / _vecSubDomains.at(0).arrBoxLength.at(1) << ","
																	 " "
		 << bbLength[2] / _vecSubDomains.at(0).arrBoxLength.at(2) << std::endl;
#endif

	std::array<double, 3> bl = _vecSubDomains.at(0).arrBoxLength;
	std::vector<Molecule>* vecParticlesVap = &_vecSubDomains.at(0).vecParticles;
	std::vector<Molecule>* vecParticlesLiq = nullptr;
	if(_nSystemType != ST_HOMOGENEOUS)
		vecParticlesLiq = &_vecSubDomains.at(1).vecParticles;

	// Count added particles
	uint64_t numAddedParticlesLocal = 0;
	uint64_t numAddedParticlesFreespaceLocal = 0;

	// set componnet
	std::vector<Component>* ptrComponents = global_simulation->getEnsemble()->getComponents();

	uint64_t bi[3];  // block index
	for(bi[0] = startIndex[0]; bi[0] < startIndex[0] + numBlocks[0]; ++bi[0]) {
		for(bi[2] = startIndex[2]; bi[2] < startIndex[2] + numBlocks[2]; ++bi[2]) {
			for(bi[1] = startIndex[1]; bi[1] < startIndex[1] + numBlocks[1]; ++bi[1]) {
				// shift particle position
				double dShift[3];
				for(uint8_t di = 0; di < 3; ++di)
					dShift[di] = bi[di] * bl.at(di);

				std::vector<Molecule>* ptrVec;
				std::vector<uint32_t>* ptrChangeVec;
				if(_nSystemType != ST_HOMOGENEOUS && bi[1] >= _nIndexLiqBeginY && bi[1] <= _nIndexLiqEndY) {
					ptrVec = vecParticlesLiq;
					ptrChangeVec = &_vecChangeCompIDsLiq;
				} else {
					ptrVec = vecParticlesVap;
					ptrChangeVec = &_vecChangeCompIDsVap;
				}

				for(auto&& mi : *ptrVec) {
					Molecule mol = mi;
					double r[3];
					for(uint8_t di = 0; di < 3; ++di) {
						r[di] = mol.r(di) + dShift[di];
						mol.setr(di, r[di]);
					}

					// Add particle to container
					double ry = r[1];
					bool bIsInsideFreespace = false;
					if(_nSystemType != ST_HOMOGENEOUS)
						bIsInsideFreespace = (ry > _fspY[0] && ry < _fspY[1]) || (ry > _fspY[2] && ry < _fspY[3]) ||
											 (ry > _fspY[4] && ry < _fspY[5]);

					if(particleContainer->isInBoundingBox(r)) {
						if(not bIsInsideFreespace) {
							// set component
							uint32_t cid = mol.componentid();
							Component* comp = &ptrComponents->at(ptrChangeVec->at(cid));
							mol.setComponent(comp);
							mol.setid(numAddedParticlesLocal);
							// inbox check already performed a few lines earlier.
							bool added = particleContainer->addParticle(mol, true, false);
							mardyn_assert(added);
							numAddedParticlesLocal++;
						} else {
							numAddedParticlesFreespaceLocal++;
						}
					}
				}
			}
		}
	}


	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(numAddedParticlesLocal);//number of local molecules
	domainDecomp->collCommScanSum();
	unsigned long idOffset = domainDecomp->collCommGetUnsLong() - numAddedParticlesLocal;
	domainDecomp->collCommFinalize();
	// fix ID's to be unique:
	for (auto mol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mol.isValid(); ++mol) {
		mol->setid(mol->getID() + idOffset);
	}


	// update global number of particles, perform number checks
	uint64_t numParticlesLocal = particleContainer->getNumberOfParticles();
	uint64_t numParticlesGlobal = 0;
	uint64_t numAddedParticlesFreespaceGlobal = 0;
	mardyn_assert(numParticlesLocal == numAddedParticlesLocal);
	domainDecomp->collCommInit(2);
	domainDecomp->collCommAppendUnsLong(numParticlesLocal);
	domainDecomp->collCommAppendUnsLong(numAddedParticlesFreespaceLocal);
	domainDecomp->collCommAllreduceSum();
	numParticlesGlobal = domainDecomp->collCommGetUnsLong();
	numAddedParticlesFreespaceGlobal = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
	mardyn_assert(numParticlesGlobal == _numParticlesTotal - numAddedParticlesFreespaceGlobal);

	Log::global_log->info() << "Number of particles calculated by number of blocks  : " << std::setw(24) << _numParticlesTotal
					   << std::endl;
	Log::global_log->info() << "Number of particles located in freespace (not added): " << std::setw(24)
					   << numAddedParticlesFreespaceGlobal << std::endl;
	Log::global_log->info() << "Number of particles added to particle container     : " << std::setw(24) << numParticlesGlobal
					   << std::endl;

	if(domainDecomp->getRank() == 0 && numParticlesGlobal != _numParticlesTotal - numAddedParticlesFreespaceGlobal) {
		Log::global_log->info() << "Number of particles: " << numParticlesGlobal << " (added)"
																			   " != "
						   << (_numParticlesTotal - numAddedParticlesFreespaceGlobal) << " (expected). Program exit ..."
						   << std::endl;
		mardyn_exit(-1);
	}

	global_simulation->timers()->stop("REPLICA_GENERATOR_VLE_INPUT");
	global_simulation->timers()->setOutputString("REPLICA_GENERATOR_VLE_INPUT", "Initial IO took:                 ");
	Log::global_log->info() << "Initial IO took:                 "
							<< global_simulation->timers()->getTime("REPLICA_GENERATOR_VLE_INPUT") << " sec"
							<< std::endl;
	Log::global_log->info() << "------------------------------------------------------------------------" << std::endl;

	return numParticlesGlobal;
}

