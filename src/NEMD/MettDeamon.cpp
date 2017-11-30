/*
 * MettDeamon.cpp
 *
 *  Created on: 03.04.2017
 *      Author: thet
 */

#include "MettDeamon.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/xmlfileUnits.h"
#include "utils/Random.h"

#include <map>
#include <array>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace std;

MettDeamon::MettDeamon()
	: 	_dDensityReservoir(0.),
		_dAreaXZ(0.),
		_dInvDensityArea(0.),
		_dY(0.),
		_dYInit(0.),
		_dYsum(0.),
		_velocityBarrier(0.),
		_dSlabWidthInit(0),
		_dSlabWidth(0),
		_dReservoirWidthY(0),
		_nUpdateFreq(0),
		_nWriteFreqRestart(0),
		_nMaxMoleculeID(0),
		_nMaxReservoirMoleculeID(0),
		_nNumMoleculesDeletedLocal(0),
		_nNumMoleculesDeletedGlobal(0),
		_nNumMoleculesDeletedGlobalAlltime(0),
		_nNumMoleculesChangedLocal(0),
		_nNumMoleculesChangedGlobal(0),
		_nNumMoleculesTooFast(0),
		_nNumMoleculesTooFastGlobal(0),
		_reservoirNumMolecules(0),
		_reservoirSlabs(0),
		_nSlabindex(0),
		_nReadReservoirMethod(RRM_UNKNOWN),
		_nMovingDirection(MD_UNKNOWN),
		_nFeedRateMethod(FRM_UNKNOWN),
		_reservoirFilename("unknown"),
		_bIsRestart(false),
		_nNumValsSummation(0),
		_numDeletedMolsSum(0),
		_dDeletedMolsPerTimestep(0.),
		_dInvNumTimestepsSummation(0.),
		_bMirrorActivated(false),
		_dMirrorPosY(0.),
		_nDeleteNonVolatile(0),
		_dMoleculeDiameter(1.0),
		_dTransitionPlanePosY(0.0),
		_dDensityTarget(0.0),
		_dVolumeCV(0.0)
{
	_dAreaXZ = global_simulation->getDomain()->getGlobalLength(0) * global_simulation->getDomain()->getGlobalLength(2);

	// init restart file
	std::ofstream ofs("MettDeamonRestart.dat", std::ios::out);
	std::stringstream outputstream;
	outputstream << "     simstep" << "   slabIndex" << "                  deltaY" << std::endl;
	ofs << outputstream.str();
	ofs.close();

	// summation of deleted molecules
	_listDeletedMolecules.clear();
	_listDeletedMolecules.push_back(0);

	// init identity change vector
	uint8_t nNumComponents = global_simulation->getEnsemble()->getComponents()->size();
	_vecChangeCompIDsFreeze.resize(nNumComponents);
	_vecChangeCompIDsUnfreeze.resize(nNumComponents);
	std::iota (std::begin(_vecChangeCompIDsFreeze), std::end(_vecChangeCompIDsFreeze), 0);
	std::iota (std::begin(_vecChangeCompIDsUnfreeze), std::end(_vecChangeCompIDsUnfreeze), 0);

	// throttle parameters
	_vecThrottleFromPosY.resize(nNumComponents);
	_vecThrottleToPosY.resize(nNumComponents);
	_vecThrottleForceY.resize(nNumComponents);
	for(uint8_t cid=0; cid<nNumComponents; ++cid)
	{
		_vecThrottleFromPosY.at(cid) = 0.;
		_vecThrottleToPosY.at(cid) = 0.;
		_vecThrottleForceY.at(cid) = 0.;
	}

	// velocity barriers
	_vecVeloctiyBarriers.resize(nNumComponents+1);

	// density values
	_vecDensityValues.clear();
}

MettDeamon::~MettDeamon()
{
}

void MettDeamon::readXML(XMLfileUnits& xmlconfig)
{
	// control
	xmlconfig.getNodeValue("control/updatefreq", _nUpdateFreq);
	xmlconfig.getNodeValue("control/writefreq", _nWriteFreqRestart);
	xmlconfig.getNodeValue("control/numvals", _nNumValsSummation);
	_dInvNumTimestepsSummation = 1. / (double)(_nNumValsSummation*_nUpdateFreq);

	// vmax
//	xmlconfig.getNodeValue("control/vmax", _velocityBarrier);
	{
		uint8_t numNodes = 0;
		XMLfile::Query query = xmlconfig.query("control/vmax");
		numNodes = query.card();
		global_log->info() << "Number of vmax values: " << (uint32_t)numNodes << endl;
		if(numNodes < 1) {
			global_log->error() << "No vmax values defined in XML-config file. Program exit ..." << endl;
			Simulation::exit(-1);
		}
		string oldpath = xmlconfig.getcurrentnodepath();
		XMLfile::Query::const_iterator nodeIter;
		for( nodeIter = query.begin(); nodeIter; nodeIter++ ) {
			xmlconfig.changecurrentnode(nodeIter);
			uint32_t cid;
			double vmax;
			xmlconfig.getNodeValue("@cid", cid);
			xmlconfig.getNodeValue(".", vmax);
			_vecVeloctiyBarriers.at(cid) = vmax;
		}
		xmlconfig.changecurrentnode(oldpath);
	}
	// Feed from left/right
	{
		_nMovingDirection = MD_UNKNOWN;
		int nVal = 0;
		xmlconfig.getNodeValue("control/direction", nVal);
		if(1 == nVal)
			_nMovingDirection = MD_LEFT_TO_RIGHT;
		else if(2 == nVal)
			_nMovingDirection = MD_RIGHT_TO_LEFT;
	}
	// Feed rate method
	{
		_nFeedRateMethod = FRM_UNKNOWN;
		int nVal = 0;
		xmlconfig.getNodeValue("control/frmethod", nVal);
		if(1 == nVal)
			_nFeedRateMethod = FRM_DELETED_MOLECULES;
		else if(2 == nVal)
			_nFeedRateMethod = FRM_CHANGED_MOLECULES;
		else if(3 == nVal)
			_nFeedRateMethod = FRM_DENSITY;
	}

	// reservoir
	bool bRet1 = xmlconfig.getNodeValue("reservoir/file", _reservoirFilename);
	bool bRet2 = xmlconfig.getNodeValue("reservoir/width", _dReservoirWidthY);
	xmlconfig.getNodeValue("reservoir/slabwidth", _dSlabWidthInit);
	if(true == bRet1 && false == bRet2)
		_nReadReservoirMethod = RRM_READ_FROM_FILE;
	else if(false == bRet1 && true == bRet2)
		_nReadReservoirMethod = RRM_READ_FROM_MEMORY;
	else
		_nReadReservoirMethod = RRM_AMBIGUOUS;

	// restart
	_bIsRestart = true;
	_bIsRestart = _bIsRestart && xmlconfig.getNodeValue("restart/slabindex", _nSlabindex);
	_bIsRestart = _bIsRestart && xmlconfig.getNodeValue("restart/deltaY", _dYsum);

	// mirror
	bool bRet = xmlconfig.getNodeValue("mirror/position", _dMirrorPosY);
	_bMirrorActivated = bRet;

	// change identity of fixed (frozen) molecules by component ID
	if(xmlconfig.changecurrentnode("changes")) {
		uint8_t numChanges = 0;
		XMLfile::Query query = xmlconfig.query("change");
		numChanges = query.card();
		global_log->info() << "Number of fixed molecules components: " << (uint32_t)numChanges << endl;
		if(numChanges < 1) {
			global_log->error() << "No component change defined in XML-config file. Program exit ..." << endl;
			Simulation::exit(-1);
		}
		string oldpath = xmlconfig.getcurrentnodepath();
		XMLfile::Query::const_iterator changeIter;
		for( changeIter = query.begin(); changeIter; changeIter++ ) {
			xmlconfig.changecurrentnode(changeIter);
			uint32_t nFrom, nTo;
			nFrom = nTo = 1;
			xmlconfig.getNodeValue("from", nFrom);
			xmlconfig.getNodeValue("to", nTo);
			_vecChangeCompIDsFreeze.at(nFrom-1) = nTo-1;
			_vecChangeCompIDsUnfreeze.at(nTo-1) = nFrom-1;
		}
		xmlconfig.changecurrentnode(oldpath);
		xmlconfig.changecurrentnode("..");
	}
	else {
		global_log->error() << "No component changes defined in XML-config file. Program exit ..." << endl;
		Simulation::exit(-1);
	}

#ifndef NDEBUG
	cout << "_vecChangeCompIDsFreeze:" << endl;
	for(uint32_t i=0; i<_vecChangeCompIDsFreeze.size(); ++i)
	{
		std::cout << i << ": " << _vecChangeCompIDsFreeze.at(i) << std::endl;
	}
	cout << "_vecChangeCompIDsUnfreeze:" << endl;
	for(uint32_t i=0; i<_vecChangeCompIDsUnfreeze.size(); ++i)
	{
		std::cout << i << ": " << _vecChangeCompIDsUnfreeze.at(i) << std::endl;
	}
#endif

	// molecule diameter
	xmlconfig.getNodeValue("diameter", _dMoleculeDiameter);

	// init _dY
	xmlconfig.getNodeValue("dyinit", _dYInit);

	// throttles
	// change identity of fixed (frozen) molecules by component ID
	if(xmlconfig.changecurrentnode("throttles")) {
		uint8_t numThrottles = 0;
		XMLfile::Query query = xmlconfig.query("throttle");
		numThrottles = query.card();
		global_log->info() << "Number of throttles: " << (uint32_t)numThrottles << endl;
		if(numThrottles < 1) {
			global_log->error() << "No throttles defined in XML-config file. Program exit ..." << endl;
			Simulation::exit(-1);
		}
		string oldpath = xmlconfig.getcurrentnodepath();
		XMLfile::Query::const_iterator throttleIter;
		for( throttleIter = query.begin(); throttleIter; throttleIter++ ) {
			xmlconfig.changecurrentnode(throttleIter);
			uint32_t cid;
			xmlconfig.getNodeValue("componentID", cid); cid--;
			xmlconfig.getNodeValue("pos/from", _vecThrottleFromPosY.at(cid));
			xmlconfig.getNodeValue("pos/to",   _vecThrottleToPosY.at(cid));
			xmlconfig.getNodeValue("force",    _vecThrottleForceY.at(cid));
			_vecThrottleForceY.at(cid) = abs(_vecThrottleForceY.at(cid)* -1.);
		}
		xmlconfig.changecurrentnode(oldpath);
		xmlconfig.changecurrentnode("..");
	}
	else {
		global_log->error() << "No throttles defined in XML-config file. Program exit ..." << endl;
		Simulation::exit(-1);
	}
}

void MettDeamon::ReadReservoir(DomainDecompBase* domainDecomp)
{
	switch(_nReadReservoirMethod)
	{
	case RRM_READ_FROM_MEMORY:
		this->ReadReservoirFromMemory(domainDecomp);
		break;
	case RRM_READ_FROM_FILE:
		this->ReadReservoirFromFile(domainDecomp);
		break;
	case RRM_UNKNOWN:
	case RRM_AMBIGUOUS:
	default:
		global_log->error() << "Unknown (or ambiguous) method to read reservoir for feature MettDeamon. Program exit ..." << endl;
		Simulation::exit(-1);
	}

	// Init reservoir slab index, unless restarting from checkpoint
	if(false == _bIsRestart)
		this->InitSlabIndex();

	// Determine max molecule id in reservoir
	this->DetermineMaxMoleculeIDs(domainDecomp);
}

void MettDeamon::ReadReservoirFromMemory(DomainDecompBase* domainDecomp)
{
	ParticleContainer* particleContainer = global_simulation->getMolecules();
	Domain* domain = global_simulation->getDomain();

	_reservoirSlabs = _dReservoirWidthY/_dSlabWidthInit;
	global_log->info() << "Mettdeamon-" << (uint32_t)_nMovingDirection << ": _reservoirSlabs=" << _reservoirSlabs << endl;
	_dSlabWidth = _dReservoirWidthY / (double)(_reservoirSlabs);
	_reservoir.resize(_reservoirSlabs);
	double Volume = _dReservoirWidthY * domain->getGlobalLength(0) * domain->getGlobalLength(2);
	this->InitTransitionPlane(global_simulation->getDomain() );
	global_log->info() << "Position of MettDeamons transition plane: " << _dTransitionPlanePosY << endl;

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM)
	{
		double y = tM->r(1);
//		if(true == this->IsBehindTransitionPlane(y) )
//			continue;

		Molecule mol(*tM);
		uint32_t nSlabindex;
		switch(_nMovingDirection)
		{
		case MD_LEFT_TO_RIGHT:
			if(y > _dReservoirWidthY)
				continue;
			nSlabindex = floor(y / _dSlabWidth);
			mol.setr(1, y - nSlabindex*_dSlabWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)
			mardyn_assert(nSlabindex < _reservoir.size() );
#ifndef NDEBUG
			cout << "y=" << y << ", mol.r(1)=" << mol.r(1) << ", >>>nSlabindex="<< nSlabindex << endl;
#endif
			_reservoir.at(nSlabindex).push_back(mol);
			break;
		case MD_RIGHT_TO_LEFT:
			if(y < (domain->getGlobalLength(1) - _dReservoirWidthY) )
				continue;
			double relPosY = y - (domain->getGlobalLength(1) - _dReservoirWidthY);
			nSlabindex = floor( relPosY / _dSlabWidth);
			cout << "y=" << y << ", relPosY=" << relPosY << ", mol.r(1)=" << mol.r(1) << ", >>>nSlabindex="<< nSlabindex << endl;
			mol.setr(1, y + (_reservoirSlabs-nSlabindex-1)*_dSlabWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)
#ifndef NDEBUG
			if(nSlabindex >= _reservoir.size() )
			{
				cout << "nSlabindex=" << nSlabindex << endl;
				cout << "y=" << y << endl;
				cout << "_reservoir.size()=" << _reservoir.size() << endl;
				cout << "_dTransitionPlanePosY=" << _dTransitionPlanePosY << endl;
			}
#endif
			mardyn_assert(nSlabindex < _reservoir.size() );
			_reservoir.at(nSlabindex).push_back(mol);
			break;
		}

		uint32_t cid = tM->componentid();
		// TODO: The following should be done by the addPartice method.
		std::vector<Component>* dcomponents = global_simulation->getEnsemble()->getComponents();
		mardyn_assert(cid < dcomponents->size() );
		dcomponents->at(cid).incNumMolecules();
	}
	std::vector<Molecule> currentReservoirSlab;
	// DEBUG
	for(uint32_t s=0; s<_reservoirSlabs; s++) {
		currentReservoirSlab = _reservoir.at(s);
		cout << "Rank["<<domainDecomp->getRank()<<"], Slab["<<s<<"]: currentReservoirSlab.size()=" << currentReservoirSlab.size() << endl;
	}
	// DEBUG
	// calc global number of particles in reservoir
	uint64_t numMoleculesLocal = 0;
	// max molecule id in reservoir (local)
	for(auto si : _reservoir)
	{
		for(auto mi : si)
			numMoleculesLocal++;
	}
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(numMoleculesLocal);
	domainDecomp->collCommAllreduceSum();
	_reservoirNumMolecules = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();

	global_log->info() << "Number of Mettdeamon Reservoirmolecules: " << _reservoirNumMolecules << endl;
	_dDensityReservoir = _reservoirNumMolecules / Volume;
	_dInvDensityArea = 1. / (_dAreaXZ * _dDensityReservoir);
	global_log->info() << "Density of Mettdeamon Reservoir: " << _dDensityReservoir << endl;
}

void MettDeamon::ReadReservoirFromFile(DomainDecompBase* domainDecomp)
{
	std::ifstream ifs;
	global_log->info() << "Opening Mettdeamon Reservoirfile " << _reservoirFilename << endl;
	ifs.open( _reservoirFilename.c_str() );
	if (!ifs.is_open()) {
		global_log->error() << "Could not open Mettdeamon Reservoirfile " << _reservoirFilename << endl;
		Simulation::exit(1);
	}
	global_log->info() << "Reading Mettdeamon Reservoirfile " << _reservoirFilename << endl;

	string token;
	vector<Component>& dcomponents = *(_simulation.getEnsemble()->getComponents());
	unsigned int numcomponents = dcomponents.size();
	string ntypestring("ICRVQD");
	enum Ndatatype { ICRVQDV, ICRVQD, IRV, ICRV } ntype = ICRVQD;

	double Xlength, Ylength, Zlength, V;
	Xlength = Ylength = Zlength = V = 1.0;
	while(ifs && (token != "NumberOfMolecules") && (token != "N"))
	{
		ifs >> token;

		if(token=="Length" || token=="L")
		{
			ifs >> Xlength >> Ylength >> Zlength;
			_reservoirSlabs = Ylength/_dSlabWidthInit;  // parameter sliceWidth
			global_log->info() << "Mettdeamon: _reservoirSlabs=" << _reservoirSlabs << endl;
			_dSlabWidth = Ylength / (double)(_reservoirSlabs);
			_reservoir.resize(_reservoirSlabs);
			V = Xlength * Ylength * Zlength;
		}
	}
	this->InitTransitionPlane(global_simulation->getDomain() );
	global_log->info() << "Position of MettDeamons transition plane: " << _dTransitionPlanePosY << endl;

	if((token != "NumberOfMolecules") && (token != "N")) {
		global_log->error() << "Expected the token 'NumberOfMolecules (N)' instead of '" << token << "'" << endl;
		Simulation::exit(1);
	}
	ifs >> _reservoirNumMolecules;
	global_log->info() << "Number of Mettdeamon Reservoirmolecules: " << _reservoirNumMolecules << endl;
	_dDensityReservoir = _reservoirNumMolecules / V;
	_dInvDensityArea = 1. / (_dAreaXZ * _dDensityReservoir);
	global_log->info() << "Density of Mettdeamon Reservoir: " << _dDensityReservoir << endl;

	streampos spos = ifs.tellg();
	ifs >> token;
	if((token=="MoleculeFormat") || (token == "M"))
	{
		ifs >> ntypestring;
		ntypestring.erase( ntypestring.find_last_not_of( " \t\n") + 1 );
		ntypestring.erase( 0, ntypestring.find_first_not_of( " \t\n" ) );

		if (ntypestring == "ICRVQDV") ntype = ICRVQDV;
		else if (ntypestring == "ICRVQD") ntype = ICRVQD;
		else if (ntypestring == "ICRV") ntype = ICRV;
		else if (ntypestring == "IRV")  ntype = IRV;
		else {
			global_log->error() << "Unknown molecule format '" << ntypestring << "'" << endl;
			Simulation::exit(1);
		}
	} else {
		ifs.seekg(spos);
	}
	global_log->info() << " molecule format: " << ntypestring << endl;

	if( numcomponents < 1 ) {
		global_log->warning() << "No components defined! Setting up single one-centered LJ" << endl;
		numcomponents = 1;
		dcomponents.resize( numcomponents );
		dcomponents[0].setID(0);
		dcomponents[0].addLJcenter(0., 0., 0., 1., 1., 1., 6., false);
	}
	double x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz, Vix, Viy, Viz;
	unsigned long id;
	unsigned int componentid;

	x=y=z=vx=vy=vz=q1=q2=q3=Dx=Dy=Dz=Vix=Viy=Viz=0.;
	q0=1.;

	for( unsigned long i = 0; i < _reservoirNumMolecules; i++ )
	{
		switch ( ntype ) {
			case ICRVQDV:
				ifs >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
					>> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz >> Vix >> Viy >> Viz;
				break;
			case ICRVQD:
				ifs >> id >> componentid >> x >> y >> z >> vx >> vy >> vz
					>> q0 >> q1 >> q2 >> q3 >> Dx >> Dy >> Dz;
				break;
			case ICRV :
				ifs >> id >> componentid >> x >> y >> z >> vx >> vy >> vz;
				break;
			case IRV :
				ifs >> id >> x >> y >> z >> vx >> vy >> vz;
				break;
		}

		if( componentid > numcomponents ) {
			global_log->error() << "Molecule id " << id << " has wrong componentid: " << componentid << ">" << numcomponents << endl;
			Simulation::exit(1);
		}
		componentid --; // TODO: Component IDs start with 0 in the program.
		//Simon changed id to i+1
		Molecule m1 = Molecule(i+1,&dcomponents[componentid],x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz);

		uint32_t nSlabindex = floor(y / _dSlabWidth);
		m1.setr(1, y - nSlabindex*_dSlabWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)

		double bbMin[3];
		double bbMax[3];
		bool bIsInsideSubdomain = false;
		domainDecomp->getBoundingBoxMinMax(global_simulation->getDomain(), bbMin, bbMax);
		bIsInsideSubdomain = x > bbMin[0] && x < bbMax[0] && y > bbMin[1] && y < bbMax[1] && z > bbMin[2] && z < bbMax[2];

		if(true == bIsInsideSubdomain)
			_reservoir.at(nSlabindex).push_back(m1);

		componentid = m1.componentid();
		// TODO: The following should be done by the addPartice method.
		dcomponents.at(componentid).incNumMolecules();

		// Print status message
		unsigned long iph = _reservoirNumMolecules / 100;
		if( iph != 0 && (i % iph) == 0 )
			global_log->info() << "Finished reading molecules: " << i/iph << "%\r" << flush;
	}

	global_log->info() << "Finished reading Mettdeamon Rerservoirmolecules: 100%" << endl;

	ifs.close();
}

void MettDeamon::DetermineMaxMoleculeIDs(DomainDecompBase* domainDecomp)
{
	ParticleContainer* particleContainer = global_simulation->getMolecules();
	uint64_t nMaxMoleculeID_local = 0;
	uint64_t nMaxReservoirMoleculeID_local = 0;

	// max molecule id in particle container (local)
	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM)
	{
		uint64_t id = tM->id();
		if(id > nMaxMoleculeID_local)
			nMaxMoleculeID_local = id;
	}

	// max molecule id in reservoir (local)
	for(auto si : _reservoir)
	{
		for(auto mi : si)
		{
			uint64_t id = mi.id();
			if(id > nMaxReservoirMoleculeID_local)
				nMaxReservoirMoleculeID_local = id;
		}
	}

	// global max IDs
#ifdef ENABLE_MPI

	MPI_Allreduce( &nMaxMoleculeID_local,          &_nMaxMoleculeID,          1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce( &nMaxReservoirMoleculeID_local, &_nMaxReservoirMoleculeID, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

#else
	_nMaxMoleculeID = nMaxMoleculeID_local;
	_nMaxReservoirMoleculeID = nMaxReservoirMoleculeID_local;
#endif
}

uint64_t MettDeamon::getnNumMoleculesDeleted2( DomainDecompBase* domainDecomposition)
{
	domainDecomposition->collCommInit(1);
		domainDecomposition->collCommAppendUnsLong(_nNumMoleculesTooFast);
		domainDecomposition->collCommAllreduceSum();
		_nNumMoleculesTooFastGlobal = domainDecomposition->collCommGetUnsLong();
		domainDecomposition->collCommFinalize();
//
//		std::cout << "Particles deleted: "<< _nNumMoleculesDeletedGlobalAlltime << std::endl;
//		std::cout << "Of which were too fast: " << _nNumMoleculesTooFastGlobal << std::endl;
		return _nNumMoleculesTooFastGlobal;
}
void MettDeamon::prepare_start(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer, double cutoffRadius)
{
	this->ReadReservoir(domainDecomp);
	double dMoleculeRadius = _dMoleculeDiameter*0.5;

	//ParticleContainer* _moleculeContainer;
	particleContainer->deleteOuterParticles();
	// fixed components

	if(true == _bIsRestart)
		return;

	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM)
	{
		double dPosY = tM->r(1);
		bool IsBehindTransitionPlane = this->IsBehindTransitionPlane(dPosY);
		if(false == IsBehindTransitionPlane)
		{
			uint32_t cid = tM->componentid();
			if(cid != _vecChangeCompIDsFreeze.at(cid))
			{
				Component* compNew = &(ptrComps->at(_vecChangeCompIDsFreeze.at(cid) ) );
				tM->setComponent(compNew);
//				cout << "cid(new) = " << tM->componentid() << endl;
			}
		}
/*
		else
//		else if(dPosY < (_dTransitionPlanePosY+dMoleculeRadius) )
		{
			particleContainer->deleteMolecule(tM->id(), tM->r(0), tM->r(1),tM->r(2), false);
//			cout << "delete: dY = " << dPosY << endl;
			particleContainer->update();
			tM  = particleContainer->iteratorBegin();
			this->IncrementDeletedMoleculesLocal();
		}

		else if(dPosY < (_dTransitionPlanePosY+_dMoleculeDiameter) )
		{
			tM->setv(1, tM->r(1)-1.);
		}
		*/
	}
	particleContainer->update();
	particleContainer->updateMoleculeCaches();
}
void MettDeamon::init_positionMap(ParticleContainer* particleContainer)
{
	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM) {

		uint64_t mid = tM->id();
		uint32_t cid = tM->componentid();

		if(cid != _vecChangeCompIDsUnfreeze.at(cid))
		{
			//savevelo
			std::array<double,10> pos;
			pos.at(0) = tM->r(0);
			pos.at(1) = tM->r(1);
			pos.at(2) = tM->r(2);
			pos.at(3) = tM->v(0);
			pos.at(4) = tM->v(1);
			pos.at(5) = tM->v(2);
			Quaternion q = tM->q();
			pos.at(6) = q.qw();
			pos.at(7) = q.qx();
			pos.at(8) = q.qy();
			pos.at(9) = q.qz();
//			_storePosition.insert ( std::pair<unsigned long, std::array<double, 3> >(mid, pos) );
			_storePosition[tM->id()] = pos;
		}
	}
}
void MettDeamon::preForce_action(ParticleContainer* particleContainer, double cutoffRadius)
{
	double dBoxY = global_simulation->getDomain()->getGlobalLength(1);
	double dMoleculeRadius = 0.5;

	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();

	Random rnd;
	double T = 80;
	double v[3];
	double v2 = 0.;
//	_dY = _dYInit;

	particleContainer->updateMoleculeCaches();

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM)
	{
		uint8_t cid = tM->componentid();
		double dY = tM->r(1);
		bool bIsTrappedMolecule = this->IsTrappedMolecule(cid);
		bool IsBehindTransitionPlane = this->IsBehindTransitionPlane(dY);
/*
		double m = tM->mass();
		double vym = sqrt(T/m);
*/
		double m = tM->mass();
		double vm2 = 3*T/m;

		v[0] = rnd.rnd();
		v[1] = rnd.rnd();
		v[2] = rnd.rnd();
		v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
//		global_log->info() << "rnd: vx=" << v[0] << ", vy=" << v[1] << ", vz=" << v[2] << ", v2=" << v2 << endl;
		double f = sqrt(vm2/v2);

		for(uint8_t dim=0; dim<3; ++dim)
			v[dim] *= f;
		v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
//		global_log->info() << "scaled: vx=" << v[0] << ", vy=" << v[1] << ", vz=" << v[2] << ", v2=" << v2 << endl;

		if(true == bIsTrappedMolecule && true == IsBehindTransitionPlane)
		{
			Component* compNew = &(ptrComps->at(_vecChangeCompIDsUnfreeze.at(cid) ) );
			tM->setComponent(compNew);
//			tM->setv(1, abs(tM->v(1) ) );
			tM->setv(0, 0.0);
			tM->setv(1, 3*_dY);
			tM->setv(2, 0.0);
			// delete fraction of non-volatile component
			/*
			if(tM->id()+1 == 1 && _nDeleteNonVolatile%10 != 0)
			{
				particleContainer->deleteMolecule(tM->id(), tM->r(0), tM->r(1),tM->r(2), false);
	//			cout << "delete: dY = " << dPosY << endl;
				particleContainer->update();
				tM  = particleContainer->iteratorBegin();
				this->IncrementDeletedMoleculesLocal();
				_nDeleteNonVolatile++;
			}
			*/
		}
		/*
		else if(dY < (_dTransitionPlanePosY+dMoleculeRadius) && dY > _dTransitionPlanePosY && tM->v(1) < 0.0 && false == bIsFrozenMolecule)
		{
			particleContainer->deleteMolecule(tM->id(), tM->r(0), tM->r(1),tM->r(2), false);
//			cout << "delete: dY = " << dPosY << endl;
			particleContainer->update();
			tM  = particleContainer->iteratorBegin();
			this->IncrementDeletedMoleculesLocal();
		}
		*/

		// reset position of fixed molecules
		std::map<unsigned long, std::array<double,10> >::iterator it;
		if(true == bIsTrappedMolecule)
		{
			it = _storePosition.find(tM->id() );
			if(it != _storePosition.end() )
			{
				if(MD_LEFT_TO_RIGHT == _nMovingDirection)
					tM->setr(1, it->second.at(1) + _dY);
				else if(MD_RIGHT_TO_LEFT == _nMovingDirection)
					tM->setr(1, it->second.at(1) - _dY);

				if(dY < _dSlabWidth)
				{
					tM->setr(0, it->second.at(0) );
					tM->setr(2, it->second.at(2) );

					// reset quaternion (orientation)
					Quaternion q(it->second.at(6),
							it->second.at(7),
							it->second.at(8),
							it->second.at(9) );
					tM->setq(q);
				}
			}
			// limit velocity of trapped molecules
//			tM->setv(0, 0.);
			tM->setv(1, 0.);
//			tM->setv(2, 0.);

			double T = 80.;
			double m = tM->mass();
			double vm2 = T/m*4/9.;
			double v2 = tM->v2();

			if(v2 > vm2)
			{
				double f = sqrt(vm2/v2);
				tM->scale_v(f);
			}
		}
		else
		{
			double v2 = tM->v2();
			double v2max = _vecVeloctiyBarriers.at(cid+1)*_vecVeloctiyBarriers.at(cid+1);

			if(v2 > v2max)
			{
				uint64_t id = tM->id();
				double dY = tM->r(1);

//				cout << "Velocity barrier for cid+1=" << cid+1 << ": " << _vecVeloctiyBarriers.at(cid+1) << endl;
//				cout << "id=" << id << ", dY=" << dY << ", v=" << sqrt(tM->v2() ) << endl;

	//			particleContainer->deleteMolecule(tM->id(), tM->r(0), tM->r(1),tM->r(2), true);
	//			_nNumMoleculesDeletedLocal++;
	//			_nNumMoleculesTooFast++;

				double f = sqrt(v2max/v2);
				tM->scale_v(f);
			}
		}

		// mirror molecules back that are on the way to pass fixed molecule region
		if(dY <= _dSlabWidth)
			tM->setv(1, abs(tM->v(1) ) );

	}  // loop over molecules

	// DEBUG
//	_dY = 0.01;
	// DEBUG
	_dYsum += _dY;

	if (_dYsum >= _dSlabWidth)
	{
		global_log->info() << "Mett-" << (uint32_t)_nMovingDirection << ": _dYsum=" << _dYsum << ", _dSlabWidth=" << _dSlabWidth << endl;
		global_log->info() << "_dSlabWidth=" << _dSlabWidth << endl;
		global_log->info() << "_dYsum=" << _dYsum << endl;
		global_log->info() << "_nSlabindex=" << _nSlabindex << endl;

		// insert actual reservoir slab / activate next reservoir slab
		this->InsertReservoirSlab(particleContainer);
	}
	particleContainer->update();
	particleContainer->updateMoleculeCaches();
}
void MettDeamon::postForce_action(ParticleContainer* particleContainer, DomainDecompBase* domainDecomposition)
{
	unsigned long nNumMoleculesLocal = 0;
	unsigned long nNumMoleculesGlobal = 0;
	double T = 80;

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM) {

		uint8_t cid = tM->componentid();
		bool bIsTrappedMolecule = this->IsTrappedMolecule(cid);
		double v2 = tM->v2();
		double dY = tM->r(1);

		if(true == bIsTrappedMolecule)
		{
			// limit velocity of trapped molecules
//			tM->setv(0, 0.);
			tM->setv(1, 0.);
//			tM->setv(2, 0.);

			double T = 80.;
			double m = tM->mass();
			double vm2 = T/m*4/9.;
			double v2 = tM->v2();

			if(v2 > vm2)
			{
				double f = sqrt(vm2/v2);
				tM->scale_v(f);
			}
		}
		else
		{
			double v2 = tM->v2();
			double v2max = _vecVeloctiyBarriers.at(cid+1)*_vecVeloctiyBarriers.at(cid+1);

			if(v2 > v2max)
			{
				uint64_t id = tM->id();
				double dY = tM->r(1);

//				cout << "Velocity barrier for cid+1=" << cid+1 << ": " << _vecVeloctiyBarriers.at(cid+1) << endl;
//				cout << "id=" << id << ", dY=" << dY << ", v=" << sqrt(tM->v2() ) << endl;

				double f = sqrt(v2max/v2);
				tM->scale_v(f);
			}
		}
		/*
		else if(v2 > _velocityBarrier*_velocityBarrier) // v2_limit)
		{
			uint64_t id = tM->id();
			double dY = tM->r(1);

			cout << "cid+1=" << cid+1 << endl;
			cout << "id=" << id << endl;
			cout << "dY=" << dY << endl;
			cout << "v2=" << v2 << endl;

			particleContainer->deleteMolecule(tM->id(), tM->r(0), tM->r(1),tM->r(2), true);
			_nNumMoleculesDeletedLocal++;
			_nNumMoleculesTooFast++;
			continue;
			this->IncrementDeletedMoleculesLocal();
		}
		*/
/*
 * limit force, velocities => repair dynamic
 *
		if(false == bIsFrozenMolecule &&
				dY > (_dTransitionPlanePosY + _vecThrottleFromPosY.at(cid) ) &&
				dY < (_dTransitionPlanePosY +   _vecThrottleToPosY.at(cid) ) )
		{
			double m = tM->mass();
//			global_log->info() << "m=" << m << endl;
			double vm = sqrt(T/m);
//			global_log->info() << "vym" << vym << endl;
			double v[3];
			for(uint8_t dim=0; dim<3; ++dim)
			{
				v[dim] = tM->v(dim);

				if(abs(v[dim]) > vm)
				{
					if(v[dim]>0.)
						v[dim] += abs(v[dim])/_vecThrottleForceY.at(cid)*(-1.);
					else
						v[dim] += abs(v[dim])/_vecThrottleForceY.at(cid);

					if(abs(v[dim]) > 5*vm)
					{
						if(v[dim]>0.)
							v[dim] = vm;
						else
							v[dim] = vm*-1;
					}
					tM->setv(dim,v[dim]);
				}
			}
//			global_log->info() << "_dTransitionPlanePosY=" << _dTransitionPlanePosY << endl;
//			global_log->info() << "dY=" << dY << endl;
		}
*/
		// mirror, to simulate VLE
		if(true == _bMirrorActivated)
		{
			if(tM->r(1) >= _dMirrorPosY)
				tM->setv(1, -1.*abs(tM->v(1) ) );
		}
	}
//	particleContainer->update();
//	particleContainer->updateMoleculeCaches();
	nNumMoleculesLocal = particleContainer->getNumberOfParticles();

	// delta y berechnen: alle x Zeitschritte
	if(global_simulation->getSimulationStep() % _nUpdateFreq == 0)
	{
		// update global number of particles / calc global number of deleted particles
		domainDecomposition->collCommInit(3);
		domainDecomposition->collCommAppendUnsLong(nNumMoleculesLocal);
		domainDecomposition->collCommAppendUnsLong(_nNumMoleculesDeletedLocal);
		domainDecomposition->collCommAppendUnsLong(_nNumMoleculesChangedLocal);
		domainDecomposition->collCommAllreduceSum();
		nNumMoleculesGlobal = domainDecomposition->collCommGetUnsLong();
		_nNumMoleculesDeletedGlobal = domainDecomposition->collCommGetUnsLong();
		_nNumMoleculesChangedGlobal = domainDecomposition->collCommGetUnsLong();
		domainDecomposition->collCommFinalize();
		_nNumMoleculesDeletedGlobalAlltime += _nNumMoleculesDeletedGlobal;
		_nNumMoleculesDeletedLocal = 0;
		_nNumMoleculesChangedLocal = 0;

		// update sum and summation list
		uint64_t numMolsDeletedOrChanged;
		if(FRM_DELETED_MOLECULES == _nFeedRateMethod)
			numMolsDeletedOrChanged = _nNumMoleculesDeletedGlobal;
		else if(FRM_CHANGED_MOLECULES == _nFeedRateMethod)
			numMolsDeletedOrChanged = _nNumMoleculesChangedGlobal;
		_numDeletedMolsSum += numMolsDeletedOrChanged;
		_numDeletedMolsSum -= _listDeletedMolecules.front();
		_listDeletedMolecules.push_back(numMolsDeletedOrChanged);
		if(_listDeletedMolecules.size() > _nNumValsSummation)
			_listDeletedMolecules.pop_front();
		else
		{
			_numDeletedMolsSum = 0;
			for(auto&& vi : _listDeletedMolecules)
				_numDeletedMolsSum += vi;
		}
		_dDeletedMolsPerTimestep = _numDeletedMolsSum * _dInvNumTimestepsSummation;
		if(FRM_DELETED_MOLECULES == _nFeedRateMethod || FRM_CHANGED_MOLECULES == _nFeedRateMethod)
			this->calcDeltaY();
		else if(FRM_DENSITY == _nFeedRateMethod)
			this->calcDeltaYbyDensity();
		global_log->info() << "_nNumMoleculesDeletedGlobal = " << _nNumMoleculesDeletedGlobal << endl;
		global_log->info() << "_nNumMoleculesChangedGlobal = " << _nNumMoleculesChangedGlobal << endl;
		global_log->info() << "_numDeletedMolsSum = " << _numDeletedMolsSum << endl;
		global_log->info() << "_dDeletedMolsPerTimestep = " << _dDeletedMolsPerTimestep << endl;
		global_log->info() << "_dY = " << _dY << endl;
	}
	else
	{
		// update global number of particles
		domainDecomposition->collCommInit(1);
		domainDecomposition->collCommAppendUnsLong(nNumMoleculesLocal);
		domainDecomposition->collCommAllreduceSum();
		nNumMoleculesGlobal = domainDecomposition->collCommGetUnsLong();
		domainDecomposition->collCommFinalize();
	}
	global_simulation->getDomain()->setglobalNumMolecules(nNumMoleculesGlobal);

	// write restart file
	this->writeRestartfile();
}

void MettDeamon::writeRestartfile()
{
	if(0 != global_simulation->getSimulationStep() % _nWriteFreqRestart)
		return;

	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();

	if(domainDecomp.getRank() != 0)
		return;

	std::ofstream ofs("MettDeamonRestart.dat", std::ios::app);
	std::stringstream outputstream;

	outputstream << setw(12) << global_simulation->getSimulationStep() << setw(12) << _nSlabindex;
	outputstream << FORMAT_SCI_MAX_DIGITS << _dYsum << std::endl;

	ofs << outputstream.str();
	ofs.close();
}

void MettDeamon::NextReservoirSlab()
{
	switch(_nMovingDirection)
	{
	case MD_LEFT_TO_RIGHT:
		_nSlabindex--;
		if(_nSlabindex == -1)
		{
			_nSlabindex = _reservoirSlabs-1;
			_nMaxMoleculeID += _nMaxReservoirMoleculeID;
		}
		break;
	case MD_RIGHT_TO_LEFT:
		_nSlabindex++;
		if( (uint32_t)(_nSlabindex) == _reservoirSlabs)
		{
			_nSlabindex = 0;
			_nMaxMoleculeID += _nMaxReservoirMoleculeID;
		}
		break;
	}
}

void MettDeamon::calcDeltaYbyDensity()
{
	double dDensityMean = 0.;
	uint32_t numVals = 0;
	for(auto&& dVal : _vecDensityValues)
	{
		dDensityMean += dVal;
		numVals++;
	}
	double dInvNumVals = 1./((double)(numVals));
	dDensityMean *= dInvNumVals;
	double dDensityDelta = _dDensityTarget - dDensityMean;
	if(dDensityDelta <= 0.)
		_dY = 0.;
	else
		_dY = dDensityDelta/_dDensityReservoir*dInvNumVals*_dVolumeCV/_dAreaXZ;
}

void MettDeamon::InsertReservoirSlab(ParticleContainer* particleContainer)
{
	cout << "INSERT: Mett-" << (uint32_t)_nMovingDirection << endl;
	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
	std::vector<Molecule> currentReservoirSlab = _reservoir.at(_nSlabindex);
	global_log->info() << "currentReservoirSlab.size()=" << currentReservoirSlab.size() << endl;

	for(auto mi : currentReservoirSlab)
	{
		uint64_t id = mi.id();
		uint32_t cid = mi.componentid();
		Component* compNew;
		if(_nMovingDirection==1)
			compNew = &(ptrComps->at(_vecChangeCompIDsFreeze.at(cid) ) );
		else
			compNew = &(ptrComps->at(3));

		mi.setid(_nMaxMoleculeID + id);
		mi.setComponent(compNew);
		mi.setr(1, mi.r(1) + _dYsum - _dSlabWidth);
		cout << "INSERT POS: " << mi.r(1) << ", cid=" << _vecChangeCompIDsFreeze.at(cid) << endl;
		particleContainer->addParticle(mi);
	}
	_dYsum -= _dSlabWidth;  // reset feed sum
	this->NextReservoirSlab();
}


