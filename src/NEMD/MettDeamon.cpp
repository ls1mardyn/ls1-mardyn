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

MettDeamon::MettDeamon(double cutoffRadius)
	: 	_rho_l(0.),
		_dAreaXZ(0.),
		_dInvDensityArea(0.),
		_dY(0.),
		_dYInit(0.),
		_dYsum(0.),
		_velocityBarrier(0.),
		_dSlabWidthInit(0),
		_dSlabWidth(0),
		_cutoffRadius(cutoffRadius),
		_nUpdateFreq(0),
		_nWriteFreqRestart(0),
		_maxId(0),
		_nNumMoleculesDeletedLocal(0),
		_nNumMoleculesDeletedGlobal(0),
		_nNumMoleculesDeletedGlobalAlltime(0),
		_nNumMoleculesTooFast(0),
		_nNumMoleculesTooFastGlobal(0),
		_reservoirNumMolecules(0),
		_reservoirSlabs(0),
		_nSlabindex(0),
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
		_dTransitionPlanePosY(0.0)
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
}

MettDeamon::~MettDeamon()
{
}

void MettDeamon::readXML(XMLfileUnits& xmlconfig)
{
	// control
	xmlconfig.getNodeValue("control/updatefreq", _nUpdateFreq);
	xmlconfig.getNodeValue("control/writefreq", _nWriteFreqRestart);
	xmlconfig.getNodeValue("control/vmax", _velocityBarrier);
	xmlconfig.getNodeValue("control/numvals", _nNumValsSummation);
	_dInvNumTimestepsSummation = 1. / (double)(_nNumValsSummation*_nUpdateFreq);

	// reservoir
	xmlconfig.getNodeValue("reservoir/file", _reservoirFilename);
	xmlconfig.getNodeValue("reservoir/slabwidth", _dSlabWidthInit);

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
	unsigned long maxid = 0; // stores the highest molecule ID found in the phase space file
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
//				_reservoirLastSlice = Ylength-((_reservoirSlices-1)*_cutoffRadius);
			_reservoir.resize(_reservoirSlabs);
			V = Xlength * Ylength * Zlength;
		}
	}

	if((token != "NumberOfMolecules") && (token != "N")) {
		global_log->error() << "Expected the token 'NumberOfMolecules (N)' instead of '" << token << "'" << endl;
		Simulation::exit(1);
	}
	ifs >> _reservoirNumMolecules;
	global_log->info() << " number of Mettdeamon Reservoirmolecules: " << _reservoirNumMolecules << endl;
	_rho_l = _reservoirNumMolecules / V;
	_dInvDensityArea = 1. / (_dAreaXZ * _rho_l);
	global_log->info() << "Density of Mettdeamon Reservoir: " << _rho_l << endl;

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

//				long temp= m1.r(1)/_cutoffRadius;
//				_reservoirIter = _reservoir.begin();
		uint32_t nSlabindex = floor(y / _dSlabWidth);
		m1.setr(1, y - nSlabindex*_dSlabWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)

		double bbMin[3];
		double bbMax[3];
		bool bIsInsideSubdomain = false;
		domainDecomp->getBoundingBoxMinMax(global_simulation->getDomain(), bbMin, bbMax);
		bIsInsideSubdomain = x > bbMin[0] && x < bbMax[0] && y > bbMin[1] && y < bbMax[1] && z > bbMin[2] && z < bbMax[2];

		if(true == bIsInsideSubdomain)
			_reservoir.at(nSlabindex).push_back(m1);

		componentid=m1.componentid();
		// TODO: The following should be done by the addPartice method.
		dcomponents[componentid].incNumMolecules();

		if(id > maxid) maxid = id;
		_maxId = maxid;

		// Print status message
		unsigned long iph = _reservoirNumMolecules / 100;
		if( iph != 0 && (i % iph) == 0 )
			global_log->info() << "Finished reading molecules: " << i/iph << "%\r" << flush;
	}

	global_log->info() << "Finished reading Mettdeamon Rerservoirmolecules: 100%" << endl;

	ifs.close();
	if(false == _bIsRestart)
		_nSlabindex = _reservoirSlabs-1;
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
	_dTransitionPlanePosY = 2*_dSlabWidth;

	//ParticleContainer* _moleculeContainer;
	particleContainer->deleteOuterParticles();
	// fixed components

	Molecule* tM;
	double dPosY;
	double dBoxY = global_simulation->getDomain()->getGlobalLength(1);
	global_log->info() << "Position of MettDeamons transition plane: " << _dTransitionPlanePosY << endl;

	if(true == _bIsRestart)
		return;

	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM)
	{
		dPosY = tM->r(1);
		if(dPosY < (_dTransitionPlanePosY) )  // -dMoleculeRadius) )
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
			std::array<double, 6> pos;
			pos.at(0) = tM->r(0);
			pos.at(1) = tM->r(1);
			pos.at(2) = tM->r(2);
			pos.at(3) = tM->v(0);
			pos.at(4) = tM->v(1);
			pos.at(5) = tM->v(2);
//			_storePosition.insert ( std::pair<unsigned long, std::array<double, 3> >(mid, pos) );
			_storePosition[tM->id()] = pos;
		}
	}
}
void MettDeamon::preForce_action(ParticleContainer* particleContainer, double cutoffRadius)
{
	double dBoxY = global_simulation->getDomain()->getGlobalLength(1);
	double dMoleculeRadius = 0.5;
	uint32_t cid;
	std::map<unsigned long, std::array<double, 6> >::iterator it;

	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();

	Random rnd;
	double T = 0.88;
	double v[3];
	double v2 = 0.;
//	_dY = _dYInit;

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM)
	{
		cid = tM->componentid();
		double dY = tM->r(1);
		bool bIsFrozenMolecule = cid != _vecChangeCompIDsUnfreeze.at(cid);

		v2 = 0.;
		v[0] = rnd.gaussDeviate(T);
		v[1] = rnd.gaussDeviate(T);
		v[2] = rnd.gaussDeviate(T);
		for(uint8_t dim=0; dim<3; ++dim)
			v2 += v[dim]*v[dim];
//		global_log->info() << "vx=" << v[0] << ", vy=" << v[1] << ", vz=" << v[2] << ", v2=" << v2 << endl;

		if(dY > (_dTransitionPlanePosY-dMoleculeRadius) && true == bIsFrozenMolecule)
		{
			Component* compNew = &(ptrComps->at(_vecChangeCompIDsUnfreeze.at(cid) ) );
			tM->setComponent(compNew);
//			tM->setv(1, abs(tM->v(1) ) );
			tM->setv(0, v[0]);
			tM->setv(1, v[1]);
			tM->setv(2, v[2]);
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
		if(true == bIsFrozenMolecule)
		{
			it = _storePosition.find(tM->id() );
			if(it != _storePosition.end() )
			{
				tM->setr(0, it->second.at(0) );
				tM->setr(1, it->second.at(1) + this->getDeltaY() );
				tM->setr(2, it->second.at(2) );
			}
			// restore velocity
			tM->setv(0,it->second.at(3) );
			tM->setv(1,it->second.at(4) );
			tM->setv(2,it->second.at(5) );
		}

	}  // loop over molecules

	_dYsum += this->getDeltaY();

	if (_dYsum >= _dSlabWidth)
	{
		global_log->info() << "_dSlabWidth = " << _dSlabWidth << endl;
		global_log->info() << "_dYsum = " << _dYsum << endl;
		global_log->info() << "_nSlabindex = " << _nSlabindex << endl;
		std::vector<Molecule> currentReservoirSlab = _reservoir.at(_nSlabindex);

		for(auto mi : currentReservoirSlab)
		{
			unsigned long tempId = mi.id();
			uint32_t cid = mi.componentid();
			Component* compNew = &(ptrComps->at(_vecChangeCompIDsFreeze.at(cid) ) );
			mi.setid(_maxId + tempId);
			mi.setComponent(compNew);
			mi.setr(1, mi.r(1) + _dYsum - _dSlabWidth);
			particleContainer->addParticle(mi);
		}

		_dYsum -= _dSlabWidth;
		_nSlabindex--;
		if(_nSlabindex == -1)
		{
			_nSlabindex = _reservoirSlabs-1;
			_maxId = _maxId + _reservoirNumMolecules;
		}
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

		uint32_t cid = tM->componentid();
		bool bIsFrozenMolecule = cid != _vecChangeCompIDsUnfreeze.at(cid);
		double v2 = tM->v2();
		double dY = tM->r(1);

		if(true == bIsFrozenMolecule)
		{
			tM->setv(0, 0.);
			tM->setv(1, 0.);
			tM->setv(2, 0.);
		}
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
//			particleContainer->update();
			continue;
			this->IncrementDeletedMoleculesLocal();
		}

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

		// mirror, to simulate VLE
		if(true == _bMirrorActivated)
		{
			if(tM->r(1) >= _dMirrorPosY)
				tM->setv(1, -1.*tM->v(1) );
		}
	}
	particleContainer->update();
	particleContainer->updateMoleculeCaches();
	nNumMoleculesLocal = particleContainer->getNumberOfParticles();

	// delta y berechnen: alle x Zeitschritte
	if(global_simulation->getSimulationStep() % _nUpdateFreq == 0)
	{
		// update global number of particles / calc global number of deleted particles
		domainDecomposition->collCommInit(2);
		domainDecomposition->collCommAppendUnsLong(nNumMoleculesLocal);
		domainDecomposition->collCommAppendUnsLong(_nNumMoleculesDeletedLocal);
		domainDecomposition->collCommAllreduceSum();
		nNumMoleculesGlobal = domainDecomposition->collCommGetUnsLong();
		_nNumMoleculesDeletedGlobal = domainDecomposition->collCommGetUnsLong();
		domainDecomposition->collCommFinalize();
		_nNumMoleculesDeletedGlobalAlltime += _nNumMoleculesDeletedGlobal;
		_nNumMoleculesDeletedLocal = 0;

		// update sum and summation list
		_numDeletedMolsSum += _nNumMoleculesDeletedGlobal;
		_numDeletedMolsSum -= _listDeletedMolecules.front();
		_listDeletedMolecules.push_back(_nNumMoleculesDeletedGlobal);
		if(_listDeletedMolecules.size() > _nNumValsSummation)
			_listDeletedMolecules.pop_front();
		else
		{
			_numDeletedMolsSum = 0;
			for(auto&& vi : _listDeletedMolecules)
				_numDeletedMolsSum += vi;
		}
		_dDeletedMolsPerTimestep = _numDeletedMolsSum * _dInvNumTimestepsSummation;
		this->calcDeltaY();
		global_log->info() << "_nNumMoleculesDeletedGlobal = " << _nNumMoleculesDeletedGlobal << endl;
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


