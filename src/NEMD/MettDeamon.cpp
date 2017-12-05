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
#ifdef ENABLE_MPI
#include "parallel/ParticleData.h"
#include "parallel/DomainDecomposition.h"
#endif
#include "particleContainer/ParticleContainer.h"
#include "utils/xmlfileUnits.h"
#include "utils/Random.h"
#include "io/ReplicaGenerator.h"  // class MoleculeDataReader

#include <map>
#include <array>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace std;

MettDeamon::MettDeamon() :
		_reservoir(nullptr),
		_bIsRestart(false),
		_bMirrorActivated(false),
		_dAreaXZ(0.),
		_dInvDensityArea(0.),
		_dY(0.),
		_dYInit(0.),
		_dYsum(0.),
		_velocityBarrier(0.),
		_dDeletedMolsPerTimestep(0.),
		_dInvNumTimestepsSummation(0.),
		_dMirrorPosY(0.),
		_dMoleculeDiameter(1.0),
		_dTransitionPlanePosY(0.0),
		_dDensityTarget(0.0),
		_dVolumeCV(0.0),
		_dFeedRate(0.0),
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
		_nMovingDirection(MD_UNKNOWN),
		_nFeedRateMethod(FRM_UNKNOWN),
		_nZone2Method(Z2M_UNKNOWN),
		_nNumValsSummation(0),
		_numDeletedMolsSum(0),
		_nDeleteNonVolatile(0)
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

	// particle reservoir
	_reservoir = new Reservoir(this);
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
		else if(4 == nVal)
		{
			_nFeedRateMethod = FRM_CONSTANT;
			xmlconfig.getNodeValue("control/feedrate", _dFeedRate);
		}
	}

	// Zone2 method
	{
		_nZone2Method = FRM_UNKNOWN;
		int nVal = 0;
		xmlconfig.getNodeValue("control/z2method", nVal);
		if(1 == nVal)
			_nZone2Method = Z2M_RESET_ALL;
		else if(2 == nVal)
			_nZone2Method = Z2M_RESET_YPOS_ONLY;
	}

	// reservoir
	{
		string oldpath = xmlconfig.getcurrentnodepath();
		xmlconfig.changecurrentnode("reservoir");
		_reservoir->readXML(xmlconfig);
		xmlconfig.changecurrentnode(oldpath);
	}

	// restart
	_bIsRestart = true;
	_bIsRestart = _bIsRestart && xmlconfig.getNodeValue("restart/binindex", _restartInfo.nBindindex);
	_bIsRestart = _bIsRestart && xmlconfig.getNodeValue("restart/deltaY", _restartInfo.dYsum);

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

//#ifndef NDEBUG
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
//#endif

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

void MettDeamon::findMaxMoleculeID(DomainDecompBase* domainDecomp)
{
	ParticleContainer* particleContainer = global_simulation->getMolecules();
	uint64_t nMaxMoleculeID_local = 0;

	// max molecule id in particle container (local)
	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM)
	{
		uint64_t id = tM->id();
		if(id > nMaxMoleculeID_local)
			nMaxMoleculeID_local = id;
	}

	// global max IDs
#ifdef ENABLE_MPI

	MPI_Allreduce( &nMaxMoleculeID_local, &_nMaxMoleculeID, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

#else
	_nMaxMoleculeID = nMaxMoleculeID_local;
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
	_reservoir->readParticleData(domainDecomp);
	_dInvDensityArea = 1. / (_dAreaXZ * _reservoir->getDensity() );
	// Activate reservoir bin with respect to restart information
	if(true == _bIsRestart)
		this->initRestart();

	this->InitTransitionPlane(global_simulation->getDomain() );
	cout << domainDecomp->getRank() << ": Position of MettDeamons transition plane: " << _dTransitionPlanePosY << endl;

	// find max molecule ID in particle container
	this->findMaxMoleculeID(domainDecomp);

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

		bool bIsTrappedMolecule = this->IsTrappedMolecule(cid);
		if(true == bIsTrappedMolecule)
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
			if(MD_LEFT_TO_RIGHT == _nMovingDirection)
				tM->setv(1, 3*_dY);
			else if(MD_RIGHT_TO_LEFT == _nMovingDirection)
				tM->setv(1, -3*_dY);
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
				{
					// reset y-position
					tM->setr(1, it->second.at(1) + _dY);
					// reset x,z-position
					if(dY < _reservoir->getBinWidth() || Z2M_RESET_ALL == _nZone2Method)
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
				else if(MD_RIGHT_TO_LEFT == _nMovingDirection)
				{
					// reset y-position
					tM->setr(1, it->second.at(1) - _dY);
					// reset x,z-position
					if(dY > (dBoxY - _reservoir->getBinWidth() ) || Z2M_RESET_ALL == _nZone2Method)
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
		if(dY <= _reservoir->getBinWidth() )
			tM->setv(1, abs(tM->v(1) ) );

	}  // loop over molecules

	// DEBUG
	if(FRM_CONSTANT == _nFeedRateMethod)
		_dY = _dFeedRate;
	// DEBUG
	_dYsum += _dY;

	if (_dYsum >= _reservoir->getBinWidth())
	{
		global_log->info() << "Mett-" << (uint32_t)_nMovingDirection << ": _dYsum=" << _dYsum << ", _dSlabWidth=" << _reservoir->getBinWidth() << endl;
		global_log->info() << "_dSlabWidth=" << _reservoir->getBinWidth() << endl;
		global_log->info() << "_dYsum=" << _dYsum << endl;
		global_log->info() << "_reservoir->getActualBinIndex()=" << _reservoir->getActualBinIndex() << endl;

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

	outputstream << setw(12) << global_simulation->getSimulationStep() << setw(12) << _reservoir->getActualBinIndex();
	outputstream << FORMAT_SCI_MAX_DIGITS << _dYsum << std::endl;

	ofs << outputstream.str();
	ofs.close();
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
		_dY = dDensityDelta/_reservoir->getDensity()*dInvNumVals*_dVolumeCV/_dAreaXZ;
}

void MettDeamon::InitTransitionPlane(Domain* domain)
{
	double dBoxY = domain->getGlobalLength(1);
	if(MD_LEFT_TO_RIGHT == _nMovingDirection)
		_dTransitionPlanePosY = 2*_reservoir->getBinWidth();
	else
		_dTransitionPlanePosY = dBoxY - 2*_reservoir->getBinWidth();
}

void MettDeamon::InsertReservoirSlab(ParticleContainer* particleContainer)
{
	cout << "INSERT: Mett-" << (uint32_t)_nMovingDirection << endl;
	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
	std::vector<Molecule>& currentReservoirSlab = _reservoir->getParticlesActualBin();
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
		mi.setr(1, mi.r(1) + _dYsum - _reservoir->getBinWidth() );
		cout << "INSERT POS: " << mi.r(1) << ", cid=" << _vecChangeCompIDsFreeze.at(cid) << endl;
		particleContainer->addParticle(mi);
	}
	_dYsum -= _reservoir->getBinWidth();  // reset feed sum
	_reservoir->nextBin(_nMaxMoleculeID);
}

void MettDeamon::initRestart()
{
	bool bRet = _reservoir->activateBin(_restartInfo.nBindindex);
	if(false == bRet)
	{
		global_log->info() << "Failed to activate reservoir bin after restart! Program exit ... " << endl;
		Simulation::exit(-1);
	}
	_dYsum = _restartInfo.dYsum;
}



// class Reservoir
Reservoir::Reservoir(MettDeamon* parent) :
	_parent(parent),
	_moleculeDataReader(nullptr),
	_binQueue(nullptr),
	_numMoleculesRead(0),
	_numMoleculesGlobal(0),
	_nMaxMoleculeID(0),
	_nMoleculeFormat(ICRVQD),
	_nReadMethod(RRM_UNKNOWN),
	_dReadWidthY(0.0),
	_dBinWidthInit(0.0),
	_dBinWidth(0.0),
	_dDensity(0.0),
	_dVolume(0.0),
	_strFilename("unknown"),
	_strFilenameHeader("unknown")
{
	// allocate BinQueue
	_binQueue = new BinQueue();

	// init identity change vector
	uint8_t nNumComponents = global_simulation->getEnsemble()->getComponents()->size();
	_vecChangeCompIDs.resize(nNumComponents);
	std::iota (std::begin(_vecChangeCompIDs), std::end(_vecChangeCompIDs), 0);
}

void Reservoir::readXML(XMLfileUnits& xmlconfig)
{
	std::string strType = "unknown";
	bool bRet1 = xmlconfig.getNodeValue("file@type", strType);
	bool bRet2 = xmlconfig.getNodeValue("width", _dReadWidthY);
	xmlconfig.getNodeValue("binwidth", _dBinWidthInit);
	_nReadMethod = RRM_UNKNOWN;
	if(true == bRet1 && false == bRet2)
	{
		if("ASCII" == strType) {
			_nReadMethod = RRM_READ_FROM_FILE;
			xmlconfig.getNodeValue("file", _strFilename);
		}
		else if("binary" == strType) {
			_nReadMethod = RRM_READ_FROM_FILE_BINARY;
			xmlconfig.getNodeValue("file/header", _strFilenameHeader);
			xmlconfig.getNodeValue("file/data", _strFilename);
		}
		else {
			global_log->error() << "Wrong file type='" << strType << "' specified. Programm exit ..." << endl;
			Simulation::exit(-1);
		}
	}
	else if(false == bRet1 && true == bRet2)
		_nReadMethod = RRM_READ_FROM_MEMORY;
	else
		_nReadMethod = RRM_AMBIGUOUS;

	// Possibly change component IDs
	if(xmlconfig.changecurrentnode("changes")) {
		uint8_t numChanges = 0;
		XMLfile::Query query = xmlconfig.query("change");
		numChanges = query.card();
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
			_vecChangeCompIDs.at(nFrom-1) = nTo-1;
		}
		xmlconfig.changecurrentnode(oldpath);
		xmlconfig.changecurrentnode("..");
	}
}

void Reservoir::readParticleData(DomainDecompBase* domainDecomp)
{
	cout << "Reservoir::readParticleData(...)" << endl;
	switch(_nReadMethod)
	{
	case RRM_READ_FROM_MEMORY:
		this->readFromMemory(domainDecomp);
		break;
	case RRM_READ_FROM_FILE:
		this->readFromFile(domainDecomp);
		break;
	case RRM_READ_FROM_FILE_BINARY:
		this->readFromFileBinary(domainDecomp);
		break;
	case RRM_UNKNOWN:
	case RRM_AMBIGUOUS:
	default:
		global_log->error() << "Unknown (or ambiguous) method to read reservoir for feature MettDeamon. Program exit ..." << endl;
		Simulation::exit(-1);
	}

	// sort particles into bins
	this->sortParticlesToBins();

	// volume, density
	this->calcNumMoleculesGlobal(domainDecomp);
//	mardyn_assert( (_numMoleculesGlobal == _numMoleculesRead) || (RRM_READ_FROM_MEMORY == _nReadMethod) );
	_dVolume = 1.; for(uint8_t dim=0; dim<3; dim++) _dVolume *= _arrBoxLength[dim];
	_dDensity = _numMoleculesGlobal / _dVolume;
	cout << "Density of Mettdeamon Reservoir: " << _dDensity << endl;
}

void Reservoir::sortParticlesToBins()
{
	DomainDecompBase* domainDecomp = &global_simulation->domainDecomposition();
	cout << domainDecomp->getRank() << ": Reservoir::sortParticlesToBins(...)" << endl;
	uint32_t numBins = _arrBoxLength[1] / _dBinWidthInit;
	_dBinWidth = _arrBoxLength.at(1) / (double)(numBins);
	cout << domainDecomp->getRank() << ": _dBinWidthInit="<<_dBinWidthInit<<endl;
	cout << domainDecomp->getRank() << ": _numBins="<<numBins<<endl;
	std::vector< std::vector<Molecule> > binVector;
	binVector.resize(numBins);
	Domain* domain = global_simulation->getDomain();
	uint32_t nBinIndex;
	for(auto&& mol:_particleVector)
	{
		// possibly change component IDs
		this->changeComponentID(mol, mol.componentid() );
		double y = mol.r(1);
		nBinIndex = floor(y / _dBinWidth);
//		cout << domainDecomp->getRank() << ": y="<<y<<", nBinIndex="<<nBinIndex<<", _binVector.size()="<<binVector.size()<<endl;
		mardyn_assert(nBinIndex < binVector.size() );
		mol.setr(1, y - nBinIndex*_dBinWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)
		switch(_parent->getMovingDirection() )
		{
		case MD_LEFT_TO_RIGHT:
			mol.setr(1, y - nBinIndex*_dBinWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)
			break;
		case MD_RIGHT_TO_LEFT:
			mol.setr(1, y - nBinIndex*_dBinWidth + (domain->getGlobalLength(1) - _dBinWidth) );  // positions in slabs related to origin (x,y,z) == (0,0,0)
			break;
		}
		binVector.at(nBinIndex).push_back(mol);
	}

	// add bin particle vectors to bin queue
	cout << domainDecomp->getRank() << ": add bin particle vectors to bin queue ..." << endl;
	switch(_parent->getMovingDirection() )
	{
		case MD_LEFT_TO_RIGHT:
			for (auto bit = binVector.rbegin(); bit != binVector.rend(); ++bit)
			{
				cout << domainDecomp->getRank() << ": (*bit).size()=" << (*bit).size() << endl;
				_binQueue->enque(*bit);
			}
			break;

		case MD_RIGHT_TO_LEFT:
			for(auto bin:binVector)
			{
				cout << domainDecomp->getRank() << ": bin.size()=" << bin.size() << endl;
				_binQueue->enque(bin);
			}
			break;
	}
	_binQueue->connectTailToHead();
}

void Reservoir::readFromMemory(DomainDecompBase* domainDecomp)
{
	ParticleContainer* particleContainer = global_simulation->getMolecules();
	Domain* domain = global_simulation->getDomain();

	_arrBoxLength.at(0) = domain->getGlobalLength(0);
	_arrBoxLength.at(1) = _dReadWidthY;
	_arrBoxLength.at(2) = domain->getGlobalLength(2);

	for (ParticleIterator tM = particleContainer->iteratorBegin();
	tM != particleContainer->iteratorEnd(); ++tM)
	{
//		if(true == this->IsBehindTransitionPlane(y) )
//			continue;

		Molecule mol(*tM);
		double y = mol.r(1);

		switch(_parent->getMovingDirection() )
		{
		case MD_LEFT_TO_RIGHT:
			if(y > _dReadWidthY) continue;
			break;
		case MD_RIGHT_TO_LEFT:
			if(y < (domain->getGlobalLength(1) - _dReadWidthY) ) continue;
			double relPosY = y - (domain->getGlobalLength(1) - _dReadWidthY);
			mol.setr(1, relPosY);  // move to origin x,y,z = 0,0,0
			break;
		}
		_particleVector.push_back(mol);
	}
}

void Reservoir::readFromFile(DomainDecompBase* domainDecomp)
{
	Domain* domain = global_simulation->getDomain();
	std::ifstream ifs;
	global_log->info() << "Opening Mettdeamon Reservoirfile " << _strFilename << endl;
	ifs.open( _strFilename.c_str() );
	if (!ifs.is_open()) {
		global_log->error() << "Could not open Mettdeamon Reservoirfile " << _strFilename << endl;
		Simulation::exit(1);
	}
	global_log->info() << "Reading Mettdeamon Reservoirfile " << _strFilename << endl;

	string token;
	vector<Component>& dcomponents = *(_simulation.getEnsemble()->getComponents());
	unsigned int numcomponents = dcomponents.size();
	string ntypestring("ICRVQD");
	enum Ndatatype { ICRVQDV, ICRVQD, IRV, ICRV } ntype = ICRVQD;

	double Xlength, Ylength, Zlength;
	while(ifs && (token != "NumberOfMolecules") && (token != "N"))
	{
		ifs >> token;

		if(token=="Length" || token=="L")
		{
			ifs >> Xlength >> Ylength >> Zlength;
			_arrBoxLength.at(0) = domain->getGlobalLength(0);
			_arrBoxLength.at(1) = Ylength;
			_arrBoxLength.at(2) = domain->getGlobalLength(2);
		}
	}

	if((token != "NumberOfMolecules") && (token != "N")) {
		global_log->error() << "Expected the token 'NumberOfMolecules (N)' instead of '" << token << "'" << endl;
		Simulation::exit(1);
	}
	ifs >> _numMoleculesRead;

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

	for( unsigned long i=0; i<_numMoleculesRead; i++ )
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
		Molecule mol = Molecule(i+1,&dcomponents[componentid],x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz);
/*
		uint32_t nSlabindex = floor(y / _dBinWidth);
		m1.setr(1, y - nSlabindex*_dBinWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)

		double bbMin[3];
		double bbMax[3];
		bool bIsInsideSubdomain = false;
		domainDecomp->getBoundingBoxMinMax(global_simulation->getDomain(), bbMin, bbMax);
		bIsInsideSubdomain = x > bbMin[0] && x < bbMax[0] && y > bbMin[1] && y < bbMax[1] && z > bbMin[2] && z < bbMax[2];

		if(true == bIsInsideSubdomain)
			_binVector.at(nSlabindex).push_back(m1);

		componentid = m1.componentid();
		// TODO: The following should be done by the addPartice method.
		dcomponents.at(componentid).incNumMolecules();
*/
		_particleVector.push_back(mol);

		// Print status message
		unsigned long iph = _numMoleculesRead / 100;
		if( iph != 0 && (i % iph) == 0 )
			global_log->info() << "Finished reading molecules: " << i/iph << "%\r" << flush;
	}

	cout << domainDecomp->getRank() << ": Finished reading Mettdeamon Rerservoirmolecules: 100%" << endl;

	ifs.close();
}

void Reservoir::readFromFileBinaryHeader()
{
	DomainDecompBase* domainDecomp = &global_simulation->domainDecomposition();
	cout << domainDecomp->getRank() << ": Reservoir::readFromFileBinaryHeader(...)" << endl;
	XMLfileUnits inp(_strFilenameHeader);

	if(false == inp.changecurrentnode("/mardyn")) {
		global_log->error() << "Could not find root node /mardyn in XML input file." << endl;
		global_log->fatal() << "Not a valid MarDyn XML input file." << endl;
		Simulation::exit(1);
	}

	bool bInputOk = true;
	double dCurrentTime = 0.;
	double dBL[3];
	uint64_t nNumMols;
	std::string strMoleculeFormat;
	bInputOk = bInputOk && inp.changecurrentnode("headerinfo");
	bInputOk = bInputOk && inp.getNodeValue("time", dCurrentTime);
	bInputOk = bInputOk && inp.getNodeValue("length/x", dBL[0] );
	bInputOk = bInputOk && inp.getNodeValue("length/y", dBL[1] );
	bInputOk = bInputOk && inp.getNodeValue("length/z", dBL[2] );
	bInputOk = bInputOk && inp.getNodeValue("number", nNumMols);
	bInputOk = bInputOk && inp.getNodeValue("format@type", strMoleculeFormat);
	double dVolume = 1;
	for(uint8_t di=0; di<3; ++di)
	{
		this->setBoxLength(di, dBL[di] );
		dVolume *= dBL[di];
	}
	_numMoleculesRead = nNumMols;
	this->setVolume(dVolume);
	this->setDensity(nNumMols / dVolume);

	if(false == bInputOk)
	{
		global_log->error() << "Content of file: '" << _strFilenameHeader << "' corrupted! Program exit ..." << endl;
		Simulation::exit(1);
	}

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

void Reservoir::readFromFileBinary(DomainDecompBase* domainDecomp)
{
	global_log->info() << "Reservoir::readFromFileBinary(...)" << endl;
	// read header
	this->readFromFileBinaryHeader();

#ifdef ENABLE_MPI
	if(domainDecomp->getRank() == 0) {
#endif
	global_log->info() << "Opening phase space file " << _strFilename << endl;
	std::ifstream ifs;
	ifs.open(_strFilename.c_str(), ios::binary | ios::in);
	if (!ifs.is_open()) {
		global_log->error() << "Could not open phaseSpaceFile " << _strFilename << endl;
		Simulation::exit(1);
	}

	global_log->info() << "Reading phase space file " << _strFilename << endl;

	vector<Component>& components = *(_simulation.getEnsemble()->getComponents());

	// Select appropriate reader
	switch (_nMoleculeFormat) {
		case ICRVQD: _moleculeDataReader = new MoleculeDataReaderICRVQD(); break;
		case ICRV: _moleculeDataReader = new MoleculeDataReaderICRV(); break;
		case IRV: _moleculeDataReader = new MoleculeDataReaderIRV(); break;
	}

	for (uint64_t pi=0; pi<_numMoleculesRead; pi++) {
		Molecule mol;
		_moleculeDataReader->read(ifs, mol, components);
		_particleVector.push_back(mol);
	}
#ifdef ENABLE_MPI
	}
#endif

	/* distribute molecules to other MPI processes */
#ifdef ENABLE_MPI
	unsigned long num_particles = _particleVector.size();
	MPI_CHECK( MPI_Bcast(&num_particles, 1, MPI_UNSIGNED_LONG, 0, domainDecomp->getCommunicator()) );

#define PARTICLE_BUFFER_SIZE  (16*1024)
	ParticleData particle_buff[PARTICLE_BUFFER_SIZE];
	int particle_buff_pos = 0;
	MPI_Datatype mpi_Particle;
	ParticleData::getMPIType(mpi_Particle);

	if(domainDecomp->getRank() == 0) {
		for(unsigned long i = 0; i < num_particles; ++i) {
			ParticleData::MoleculeToParticleData(particle_buff[particle_buff_pos], _particleVector[i]);
			particle_buff_pos++;
			if ((particle_buff_pos >= PARTICLE_BUFFER_SIZE) || (i == num_particles - 1)) {
				global_log->debug() << "broadcasting(sending) particles" << endl;
				MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, 0, domainDecomp->getCommunicator());
				particle_buff_pos = 0;
			}
		}
	} else {
		for(unsigned long i = 0; i < num_particles; ++i) {
			if(i % PARTICLE_BUFFER_SIZE == 0) {
				global_log->debug() << "broadcasting(receiving) particles" << endl;
				MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, 0, domainDecomp->getCommunicator());
				particle_buff_pos = 0;
			}
			Molecule mol;
			ParticleData::ParticleDataToMolecule(particle_buff[particle_buff_pos], mol);
			particle_buff_pos++;
			_particleVector.push_back(mol);
		}
	}
	global_log->debug() << "broadcasting(sending/receiving) particles complete" << endl;
#endif
	cout << domainDecomp->getRank() << ": Reading Molecules done" << endl;
}

uint64_t Reservoir::calcNumMoleculesGlobal(DomainDecompBase* domainDecomp)
{
	uint64_t numMoleculesLocal = this->getNumMoleculesLocal();
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(numMoleculesLocal);
	domainDecomp->collCommAllreduceSum();
	_numMoleculesGlobal = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();

	global_log->info() << "Number of Mettdeamon Reservoirmolecules: " << _numMoleculesGlobal << endl;
	return _numMoleculesGlobal;
}

void Reservoir::changeComponentID(Molecule& mol, const uint32_t& cid)
{
	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
	Component* compNew = &(ptrComps->at(_vecChangeCompIDs.at(cid) ) );
	mol.setComponent(compNew);
}

// queue methods
uint32_t Reservoir::getActualBinIndex() {return _binQueue->getActualBinIndex();}
uint64_t Reservoir::getNumMoleculesLocal() {return _binQueue->getNumParticles();}
uint32_t Reservoir::getNumBins() {return _binQueue->getNumBins();}
std::vector<Molecule>& Reservoir::getParticlesActualBin() {return _binQueue->getParticlesActualBin();}
void Reservoir::nextBin(uint64_t& nMaxID) {_binQueue->next(); nMaxID += _numMoleculesGlobal;}
uint64_t Reservoir::getMaxMoleculeID() {return _binQueue->getMaxID();}
bool Reservoir::activateBin(uint32_t nBinIndex){return _binQueue->activateBin(nBinIndex);}
