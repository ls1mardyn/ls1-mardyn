#include "ReplicaGeneratorVLE.h"

/* the following macro has to be defined to use math constants in cmath */
#define _USE_MATH_DEFINES  1
#include <cmath>
#include <vector>
#include <cstdint>

#include "Domain.h"

#include "parallel/DomainDecompBase.h"
#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"

using Log::global_log;

enum MoleculeFormat : uint32_t {
	ICRVQD, IRV, ICRV
};

ReplicaGeneratorVLE::ReplicaGeneratorVLE()
	:
	_numParticlesLiq(0),
	_numParticlesVap(0),
	_numParticlesTotal(0),
	_numBlocksXZ(0),
	_numBlocksLiqY(0),
	_numBlocksVapY(0),
	_nIndexLiqBeginY(0),
	_nIndexLiqEndY(0),
	_strFilePathHeaderLiq("unknown"),
	_strFilePathDataLiq("unknown"),
	_strFilePathHeaderVap("unknown"),
	_strFilePathDataVap("unknown"),
	_dBoxLengthLiqXYZ(0.0),
	_dBoxLengthVapXYZ(0.0),
	_dBoxLengthXYZ(0.0),
	_nMoleculeFormat(ICRVQD),
	_moleculeDataReader(nullptr),
	_nMaxID(0),
	_dDensityLiq(0.0),
	_dBoxVolumeLiq(0.0)
{
}

ReplicaGeneratorVLE::~ReplicaGeneratorVLE()
{
}

void ReplicaGeneratorVLE::readReplicaPhaseSpaceHeader(const std::string& strFilePathHeader, uint64_t& numParticles, double& dBoxLengthXYZ)
{
	XMLfileUnits inp(strFilePathHeader);

	if(false == inp.changecurrentnode("/mardyn")) {
		global_log->error() << "Could not find root node /mardyn in XML input file." << endl;
		global_log->fatal() << "Not a valid MarDyn XML input file." << endl;
		Simulation::exit(1);
	}

	bool bInputOk = true;
	double dCurrentTime = 0.;
	double dBoxLength[3] = {0., 0., 0.};
	uint64_t numMolecules = 0;
	std::string strMoleculeFormat;
	bInputOk = bInputOk && inp.changecurrentnode("headerinfo");
	bInputOk = bInputOk && inp.getNodeValue("time", dCurrentTime);
	bInputOk = bInputOk && inp.getNodeValue("length/x", dBoxLength[0] );
	bInputOk = bInputOk && inp.getNodeValue("length/y", dBoxLength[1] );
	bInputOk = bInputOk && inp.getNodeValue("length/z", dBoxLength[2] );
	bInputOk = bInputOk && inp.getNodeValue("number", numParticles);
	bInputOk = bInputOk && inp.getNodeValue("format@type", strMoleculeFormat);

	if(false == bInputOk)
	{
		global_log->error() << "Content of file: '" << strFilePathHeader << "' corrupted! Program exit ..." << endl;
		Simulation::exit(1);
	}

	// Box length
	if(dBoxLength[0] != dBoxLength[1] || dBoxLength[0] != dBoxLength[2] || dBoxLength[1] != dBoxLength[2])
	{
		global_log->error() << "System is not a cube! Program exit ..." << endl;
		Simulation::exit(1);
	}
	dBoxLengthXYZ = dBoxLength[0];

	if("ICRVQD" == strMoleculeFormat)
		_nMoleculeFormat = ICRVQD;
	else if("IRV" == strMoleculeFormat)
		_nMoleculeFormat = IRV;
	else if("ICRV" == strMoleculeFormat)
		_nMoleculeFormat = ICRV;
	else
	{
		global_log->error() << "Not a valid molecule format: " << strMoleculeFormat << ", program exit ..." << endl;
		Simulation::exit(1);
	}
}

void ReplicaGeneratorVLE::readReplicaPhaseSpaceData(const std::string& strFilePathData, const uint64_t& numParticles, std::vector<Molecule>& vecParticles)
{
	// Open file
	std::ifstream ifs;
	global_log->info() << "Opening phase space file " << strFilePathData << endl;
	ifs.open(strFilePathData.c_str(),
			ios::binary | ios::in);
	if (!ifs.is_open()) {
		global_log->error() << "Could not open phaseSpaceFile " << strFilePathData << endl;
		Simulation::exit(1);
	}
	global_log->info() << "Reading phase space file " << strFilePathData << endl;

	vector<Component>& components = *(_simulation.getEnsemble()->getComponents());
	unsigned int numComponents = components.size();
	unsigned long maxid = 0; // stores the highest molecule ID found in the phase space file

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

	for (uint64_t pi=0; pi<numParticles; pi++)
	{
		Molecule mol;
		_moleculeDataReader->read(ifs, mol, components);
		vecParticles.push_back(mol);

		// Print status message
		unsigned long iph = numParticles / 100;
		if (iph != 0 && (pi % iph) == 0)
			global_log->info() << "Finished reading molecules: " << pi / iph
					<< "%\r" << std::flush;
	}

	global_log->info() << "Finished reading molecules: 100%" << endl;
	global_log->info() << "Reading Molecules done" << endl;
}

void ReplicaGeneratorVLE::readXML(XMLfileUnits& xmlconfig)
{
	xmlconfig.getNodeValue("files/liquid/header", _strFilePathHeaderLiq);
	xmlconfig.getNodeValue("files/liquid/data", _strFilePathDataLiq);
	xmlconfig.getNodeValue("files/vapor/header", _strFilePathHeaderVap);
	xmlconfig.getNodeValue("files/vapor/data", _strFilePathDataVap);

	xmlconfig.getNodeValue("numblocks/xz",     _numBlocksXZ);
	xmlconfig.getNodeValue("numblocks/liquid", _numBlocksLiqY);
	xmlconfig.getNodeValue("numblocks/vapor",  _numBlocksVapY);

	// liquid blocks begin/end index
	_nIndexLiqBeginY = _numBlocksVapY;
	_nIndexLiqEndY = _numBlocksVapY + _numBlocksLiqY - 1;

	xmlconfig.getNodeValue("diameter",  _dMoleculeDiameter);

	//init
	this->init(xmlconfig);
}

void ReplicaGeneratorVLE::init(XMLfileUnits& xmlconfig)
{
	DomainDecompBase* domainDecomp = &global_simulation->domainDecomposition();
	global_log->info() << domainDecomp->getRank() << ": Init Replica VLE ..." << endl;

	// Read liquid system
	this->readReplicaPhaseSpaceHeader(_strFilePathHeaderLiq, _numParticlesLiq, _dBoxLengthLiqXYZ);
	this->readReplicaPhaseSpaceData(_strFilePathDataLiq, _numParticlesLiq, _vecParticlesLiq);

	// Read vapor system
	this->readReplicaPhaseSpaceHeader(_strFilePathHeaderVap, _numParticlesVap, _dBoxLengthVapXYZ);
	this->readReplicaPhaseSpaceData(_strFilePathDataVap, _numParticlesVap, _vecParticlesVap);

	// Box length
	if(_dBoxLengthLiqXYZ != _dBoxLengthVapXYZ)
	{
		global_log->error() << "Box length of liquid and vapor system differ! Program exit ..." << endl;
		Simulation::exit(1);
	}
	_dBoxLengthXYZ = _dBoxLengthLiqXYZ;

	// total num particles, maxID
	_numParticlesTotal = (2*_numBlocksVapY*_numParticlesVap + _numBlocksLiqY*_numParticlesLiq) * _numBlocksXZ * _numBlocksXZ;
	global_simulation->getDomain()->setglobalNumMolecules(_numParticlesTotal);
	_dBoxVolumeLiq = _dBoxLengthLiqXYZ*_dBoxLengthLiqXYZ*_dBoxLengthLiqXYZ;
	_dDensityLiq = _numParticlesLiq / _dBoxVolumeLiq;

	// update global length of domain
	double dLength[3];
	dLength[1] = (2*_numBlocksVapY + _numBlocksLiqY) * _dBoxLengthXYZ;
	dLength[0] = dLength[2] = _numBlocksXZ * _dBoxLengthXYZ;
	for(uint8_t di=0; di<3; ++di)
		global_simulation->getDomain()->setGlobalLength(di, dLength[di]);
	global_log->info() << "Domain box length = " << dLength[0] << ", " << dLength[1] << ", " << dLength[2] << endl;

/*
	// Reset domain decomposition
	if (domainDecomp != nullptr) {
		delete domainDecomp;
	}
#ifndef ENABLE_MPI
	global_log->info() << "Initializing the alibi domain decomposition ... " << endl;
	domainDecomp = new DomainDecompBase();
#else
	global_log->info() << "Initializing the standard domain decomposition ... " << endl;
	domainDecomp = (DomainDecompBase*) new DomainDecomposition();
#endif
	global_log->info() << "Initialization done" << endl;
	domainDecomp->readXML(xmlconfig);
*/

	// calc free space positions
	double dPhaseLengthVapY = _numBlocksVapY * _dBoxLengthXYZ;
	double dPhaseLengthLiqY = _numBlocksLiqY * _dBoxLengthXYZ;
	_fspY[0] = dPhaseLengthVapY - _dMoleculeDiameter;
	_fspY[1] = dPhaseLengthVapY;
	_fspY[2] = dPhaseLengthVapY + dPhaseLengthLiqY;
	_fspY[3] = dPhaseLengthVapY + dPhaseLengthLiqY + _dMoleculeDiameter;
	_fspY[4] = dPhaseLengthVapY + dPhaseLengthLiqY + dPhaseLengthVapY - _dMoleculeDiameter;
	_fspY[5] = dPhaseLengthVapY + dPhaseLengthLiqY + dPhaseLengthVapY;
}

long unsigned int ReplicaGeneratorVLE::readPhaseSpace(ParticleContainer* particleContainer,
		list<ChemicalPotential>* /*lmu*/, Domain* domain, DomainDecompBase* domainDecomp)
{
//	if(3 != domainDecomp->getRank() )
//		return 0;

	global_log->info() << domainDecomp->getRank() << ": Constructing Replica VLE ..." << endl;

	double bbMin[3];
	double bbMax[3];
	double bbLength[3];
	uint64_t numBlocks[3];
	uint64_t startIndex[3];

	for(uint8_t di=0; di<3; ++di)
	{
		bbMin[di] = domainDecomp->getBoundingBoxMin(di, domain);
		bbMax[di] = domainDecomp->getBoundingBoxMax(di, domain);
		bbLength[di] = bbMax[di] - bbMin[di];
		numBlocks[di]  =  ceil(bbLength[di] / _dBoxLengthXYZ);
		startIndex[di] = floor(bbMin[di]    / _dBoxLengthXYZ);
	}

	// Init maxID
	int nRank = domainDecomp->getRank();
	double dVolumeSubdomain = bbLength[0] * bbLength[1] * bbLength[2];
	uint64_t numParticlesPerSubdomainMax = (uint64_t) ceil(_dDensityLiq * dVolumeSubdomain);
	_nMaxID = 1 + numParticlesPerSubdomainMax * nRank;

#ifndef NDEBUG
	cout << domainDecomp->getRank() << ": numParticlesPerSubdomainMax = " << numParticlesPerSubdomainMax << endl;
	cout << domainDecomp->getRank() << ": _nMaxID (init) = " << _nMaxID << endl;
	cout << domainDecomp->getRank() << ": bbMin = " << bbMin[0] << ", " << bbMin[1] << ", " << bbMin[2] << endl;
	cout << domainDecomp->getRank() << ": bbMax = " << bbMax[0] << ", " << bbMax[1] << ", " << bbMax[2] << endl;
	cout << domainDecomp->getRank() << ": bbLength = " << bbLength[0] << ", " << bbLength[1] << ", " << bbLength[2] << endl;
	cout << domainDecomp->getRank() << ": numBlocks = " << numBlocks[0] << ", " << numBlocks[1] << ", " << numBlocks[2] << endl;
	cout << domainDecomp->getRank() << ": startIndex = " << startIndex[0] << ", " << startIndex[1] << ", " << startIndex[2] << endl;
	cout << domainDecomp->getRank() << ": _dBoxLengthXYZ = " << _dBoxLengthXYZ << endl;
	cout << domainDecomp->getRank() << ": bbLength/_dBoxLengthXYZ = " << bbLength[0]/_dBoxLengthXYZ << ", " << bbLength[1]/_dBoxLengthXYZ << ", " << bbLength[2]/_dBoxLengthXYZ << endl;
#endif

	uint64_t bi[3];  // block index
	for(bi[0]=startIndex[0]; bi[0]<startIndex[0]+numBlocks[0]; ++bi[0])
	{
		for(bi[2]=startIndex[2]; bi[2]<startIndex[2]+numBlocks[2]; ++bi[2])
		{
			for(bi[1]=startIndex[1]; bi[1]<startIndex[1]+numBlocks[1]; ++bi[1])
			{
			#ifndef NDEBUG
				cout << domainDecomp->getRank() << ": " << bi[0] << ", " << bi[1] << ", " << bi[2] << endl;
			#endif

				// shift particle position
				double dShift[3];
				for(uint8_t di=0; di<3; ++di)
					dShift[di] = bi[di] * _dBoxLengthXYZ;

				std::vector<Molecule>* ptrVec;
				if(bi[1] >= _nIndexLiqBeginY && bi[1] <= _nIndexLiqEndY)
					ptrVec = &_vecParticlesLiq;
				else
					ptrVec = &_vecParticlesVap;

				for(auto&& mi : *ptrVec)
				{
					Molecule mol = mi;
					double r[3];
					for(uint8_t di=0; di<3; ++di)
					{
						r[di] = mol.r(di) + dShift[di];
						mol.setr(di, r[di]);
					}
					mol.setid(_nMaxID);

					// Add particle to container
					double ry = r[1];
					bool bIsInsideFreespace = (ry > _fspY[0] && ry < _fspY[1]) || (ry > _fspY[2] && ry < _fspY[3]) || (ry > _fspY[4] && ry < _fspY[5]);

					if(true == particleContainer->isInBoundingBox(r) && false == bIsInsideFreespace)
					{
						particleContainer->addParticle(mol);
						_nMaxID++;
					}
				}
			}
		}
	}

	// update global number of particles
	uint64_t numParticlesLocal = particleContainer->getNumberOfParticles();
	uint64_t numParticlesGlobal = 0;
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(numParticlesLocal);
	domainDecomp->collCommAllreduceSum();
	numParticlesGlobal = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
	domain->setglobalNumMolecules(numParticlesGlobal);

	return 0;
}

