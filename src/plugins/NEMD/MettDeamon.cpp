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
#include "utils/FileUtils.h"
#include "io/ReplicaGenerator.h"  // class MoleculeDataReader

#include <map>
#include <array>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <cstdint>
#include <random>
#include <ctime>

using namespace std;

template < typename T > void shuffle( std::list<T>& lst ) // shuffle contents of a list
{
    // create a vector of (wrapped) references to elements in the list
    // http://en.cppreference.com/w/cpp/utility/functional/reference_wrapper
    std::vector< std::reference_wrapper< const T > > vec( lst.begin(), lst.end() ) ;

    // shuffle (the references in) the vector
    std::shuffle( vec.begin(), vec.end(), std::mt19937{ std::random_device{}() } ) ;

    // copy the shuffled sequence into a new list
    std::list<T> shuffled_list {  vec.begin(), vec.end() } ;

    // swap the old list with the shuffled list
    lst.swap(shuffled_list) ;
}

void create_rand_vec_ones(const uint64_t& nCount, const double& percent, std::vector<int>& v)
{
	v.resize(nCount);
	if(nCount*percent < 1)
		return;
	std::fill (v.begin(),v.end(),0);
	std::fill (v.begin(),v.begin()+static_cast<int>(round(nCount*percent)),1);

	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(v.begin(), v.end(), g);
}

void update_velocity_vectors(Random* rnd, const uint64_t& numSamples, const double&T, const double&D, const double&v_neg, const double&e_neg,
		std::vector<double>& vxi, std::vector<double>& vyi, std::vector<double>& vzi)
{
	double stdDev = sqrt(T);
	vector<double> v2xi(numSamples);
	vector<double> v2yi(numSamples);
	vector<double> v2zi(numSamples);
	vxi.resize(numSamples);
	vyi.resize(numSamples);
	vzi.resize(numSamples);
	std::fill (vxi.begin(),vxi.end(),0);
	std::fill (vyi.begin(),vyi.end(),0);
	std::fill (vzi.begin(),vzi.end(),0);

	for(uint64_t i=0; i<numSamples; i++)
	{
		double vx = rnd->gaussDeviate(stdDev);
		double vy = 1.;
		while(vy > 0.)
			vy = rnd->gaussDeviate(stdDev) + D;
		double vz = rnd->gaussDeviate(stdDev);

		// store velocities
		vxi.at(i) = vx;
		vyi.at(i) = vy;
		vzi.at(i) = vz;

		// store squared velocities
		v2xi.at(i) = vx*vx;
		v2yi.at(i) = vy*vy;
		v2zi.at(i) = vz*vz;
	}

	// --> DRIFT
	// calc drift
	double sum_vxi = std::accumulate(vxi.begin(), vxi.end(), 0.0);
	double sum_vyi = std::accumulate(vyi.begin(), vyi.end(), 0.0);
	double sum_vzi = std::accumulate(vzi.begin(), vzi.end(), 0.0);

	double dInvNumSamples = 1. / static_cast<double>(numSamples);
	double vxd = sum_vxi * dInvNumSamples;
	double vyd = sum_vyi * dInvNumSamples;
	double vzd = sum_vzi * dInvNumSamples;

	// correct drift
	for(double & it : vxi)
		it += (0.0 - vxd);

	for(double & it : vyi)
		it += (v_neg - vyd);

	for(double & it : vzi)
		it += (0.0 - vzd);

	// calc drift again
	sum_vxi = std::accumulate(vxi.begin(), vxi.end(), 0.0);
	sum_vyi = std::accumulate(vyi.begin(), vyi.end(), 0.0);
	sum_vzi = std::accumulate(vzi.begin(), vzi.end(), 0.0);

	vxd = sum_vxi * dInvNumSamples;
	vyd = sum_vyi * dInvNumSamples;
	vzd = sum_vzi * dInvNumSamples;

	// update v2 vectors
	for(uint64_t i=0; i<numSamples; i++)
	{
		double vx = vxi.at(i);
		double vy = vyi.at(i);
		double vz = vzi.at(i);

		v2xi.at(i) = vx*vx;
		v2yi.at(i) = vy*vy;
		v2zi.at(i) = vz*vz;
	}
	// <-- DRIFT


	// --> EKIN
	// calc drift
	double sum_v2xi = std::accumulate(v2xi.begin(), v2xi.end(), 0.0);
	double sum_v2yi = std::accumulate(v2yi.begin(), v2yi.end(), 0.0);
	double sum_v2zi = std::accumulate(v2zi.begin(), v2zi.end(), 0.0);

	double v2x = sum_v2xi * dInvNumSamples;
	double v2y = sum_v2yi * dInvNumSamples;
	double v2z = sum_v2zi * dInvNumSamples;

	// correct ekin
	double scale_vx = sqrt(numSamples*T/sum_v2xi);
	for(double & it : vxi)
		it *= scale_vx;

	double scale_vy = sqrt(numSamples*e_neg/sum_v2yi);
	for(double & it : vyi)
		it *= scale_vy;

	double scale_vz = sqrt(numSamples*T/sum_v2zi);
	for(double & it : vzi)
		it *= scale_vz;

	// update v2 vectors
	for(uint64_t i=0; i<numSamples; i++)
	{
		double vx = vxi.at(i);
		double vy = vyi.at(i);
		double vz = vzi.at(i);

		v2xi.at(i) = vx*vx;
		v2yi.at(i) = vy*vy;
		v2zi.at(i) = vz*vz;
	}

	// calc ekin again
	sum_v2xi = std::accumulate(v2xi.begin(), v2xi.end(), 0.0);
	sum_v2yi = std::accumulate(v2yi.begin(), v2yi.end(), 0.0);
	sum_v2zi = std::accumulate(v2zi.begin(), v2zi.end(), 0.0);

	v2x = sum_v2xi * dInvNumSamples;
	v2y = sum_v2yi * dInvNumSamples;
	v2z = sum_v2zi * dInvNumSamples;
	// <-- EKIN

	// calc drift again
	sum_vxi = std::accumulate(vxi.begin(), vxi.end(), 0.0);
	sum_vyi = std::accumulate(vyi.begin(), vyi.end(), 0.0);
	sum_vzi = std::accumulate(vzi.begin(), vzi.end(), 0.0);

	vxd = sum_vxi * dInvNumSamples;
	vyd = sum_vyi * dInvNumSamples;
	vzd = sum_vzi * dInvNumSamples;
}

MettDeamon::MettDeamon() :
		_bIsRestart(false),
		_bInitFeedrateLog(false),
		_bInitRestartLog(false),
		_dAreaXZ(0.),
		_dInvDensityArea(0.),
		_dDeletedMolsPerTimestep(0.),
		_dInvNumTimestepsSummation(0.),
		_dTransitionPlanePosY(0.0),
		_dDensityTarget(0.0),
		_dVolumeCV(0.0),
		_nUpdateFreq(0),
		_nWriteFreqRestart(0),
		_nMaxReservoirMoleculeID(0),
		_nNumMoleculesDeletedGlobalAlltime(0),
		_nMovingDirection(MD_UNKNOWN),
		_nFeedRateMethod(FRM_UNKNOWN),
		_nZone2Method(Z2M_UNKNOWN),
		_nNumValsSummation(0),
		_numDeletedMolsSum(0),
		_nDeleteNonVolatile(0)
{
	_dAreaXZ = global_simulation->getDomain()->getGlobalLength(0) * global_simulation->getDomain()->getGlobalLength(2);

	// summation of deleted molecules
	_listDeletedMolecules.clear();
	_listDeletedMolecules.push_back(0);

	// init identity change vector
	uint8_t nNumComponents = global_simulation->getEnsemble()->getComponents()->size();
	_vecChangeCompIDsFreeze.resize(nNumComponents);
	_vecChangeCompIDsUnfreeze.resize(nNumComponents);
	std::iota (std::begin(_vecChangeCompIDsFreeze), std::end(_vecChangeCompIDsFreeze), 0);
	std::iota (std::begin(_vecChangeCompIDsUnfreeze), std::end(_vecChangeCompIDsUnfreeze), 0);

	// density values
	_vecDensityValues.clear();

	// particle reservoir
	_reservoir.reset(new Reservoir(this));

	// init manipulation count for determination of the feed rate
	_feedrate.numMolecules.inserted.local = 0;
	_feedrate.numMolecules.deleted.local = 0;
	_feedrate.numMolecules.changed_to.local = 0;
	_feedrate.numMolecules.changed_from.local = 0;

	// init feed rate sum
	_feedrate.feed.sum = 0;

	// init released counter
	_released.count.local = 0;
	_released.deleted.local = 0;
	_released.init_file = false;

	// seed rand()
	srand (static_cast <unsigned> (time(0)));

	// init released velocities vector
	_released.log_v.clear();
	_released.log_v.resize(0);
	_released.init_file_vel = false;
	_released.log_freq_vel = 5000;

	// random numbers
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	int nRank = domainDecomp.getRank();
	_rnd.reset(new Random(8624+nRank));
}

void MettDeamon::readXML(XMLfileUnits& xmlconfig)
{
	// control
	xmlconfig.getNodeValue("control/updatefreq", _nUpdateFreq);
	xmlconfig.getNodeValue("control/logfreqfeed", _feedrate.log_freq);
	xmlconfig.getNodeValue("control/logfreqreleased", _released.log_freq);
	xmlconfig.getNodeValue("control/logfreqreleased_vel", _released.log_freq_vel);
	xmlconfig.getNodeValue("control/writefreq", _nWriteFreqRestart);
	xmlconfig.getNodeValue("control/numvals", _nNumValsSummation);
	_dInvNumTimestepsSummation = 1. / (double)(_nNumValsSummation*_nUpdateFreq);

	// Feed rate
	{
		xmlconfig.getNodeValue("control/feed/targetID", _feedrate.cid_target);
		xmlconfig.getNodeValue("control/feed/init", _feedrate.feed.init);

		_nMovingDirection = MD_UNKNOWN;
		int nVal = 0;
		xmlconfig.getNodeValue("control/feed/direction", nVal);
		if(1 == nVal)
			_nMovingDirection = MD_LEFT_TO_RIGHT;
		else if(2 == nVal)
			_nMovingDirection = MD_RIGHT_TO_LEFT;

		_nFeedRateMethod = FRM_UNKNOWN;
		nVal = 0;
		xmlconfig.getNodeValue("control/feed/method", nVal);
		if(1 == nVal) {
			_nFeedRateMethod = FRM_DELETED_MOLECULES;
			global_log->info() << "[MettDeamon] Feed method 1" << std::endl;
		}
		else if(2 == nVal) {
			_nFeedRateMethod = FRM_CHANGED_MOLECULES;
			global_log->info() << "[MettDeamon] Feed method 2" << std::endl;
		}
		else if(3 == nVal) {
			_nFeedRateMethod = FRM_DENSITY;
			global_log->info() << "[MettDeamon] Feed method 3" << std::endl;
		}
		else if(4 == nVal) {
			_nFeedRateMethod = FRM_CONSTANT;
			xmlconfig.getNodeValue("control/feed/target", _feedrate.feed.init);
			global_log->info() << "[MettDeamon] Feed method 4: Using constant feed rate with v = " << _feedrate.feed.init << std::endl;
		}
		else if(5 ==nVal) {
			global_log->info() << "[MettDeamon] Feed method 5: Getting feed rate from MettDeamonFeedrateDirector" << std::endl;
		}

		_feedrate.release_velo.method = RVM_UNKNOWN;
		nVal = 0;
		xmlconfig.getNodeValue("control/feed/release_velo/method", nVal);
		if(1 == nVal)
			_feedrate.release_velo.method = RVM_UNCHANGED;
		else if(2 == nVal || 3 == nVal)
		{
			_feedrate.release_velo.method = RVM_FIX_VALUE;
			xmlconfig.getNodeValue("control/feed/release_velo/fix_value", _feedrate.release_velo.fix_value);
		}
		else if (4 == nVal)
		{
			_feedrate.release_velo.method = RVM_NORM_DISTR;
			bool bRet = true;
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/vxz", _norm.fname.vxz);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/vy",  _norm.fname.vy);
			if(bRet) {  // TODO: move this to method: init()? Has to be called before method: afterForces(), within method Simulation::prepare_start()
				global_log->info() << "[MettDeamon] Release velocities uses MB from files." << std::endl;
				this->readNormDistr();
				shuffle(_norm.vxz);  // sequence should differ between processes
				shuffle(_norm.vy);   // same here
				DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
				int nRank = domainDecomp.getRank();
			}
		}
		else if (5 == nVal)
		{
			_feedrate.release_velo.method = RVM_NORM_DISTR_GENERATOR;
			double T, D;  // T:Temperature, D:Drift
			double a_neg, a_pos, v_neg, v_pos, e_neg, e_pos;
			std::string fpath;
			uint64_t numSamples;
			bool bRet = true;
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/temperature", T);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/drift", D);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/a_neg", a_neg);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/a_pos", a_pos);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/v_neg", v_neg);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/v_pos", v_pos);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/e_neg", e_neg);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/e_pos", e_pos);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/normMB/fpath", fpath);
			bRet = bRet && xmlconfig.getNodeValue("control/feed/release_velo/numSamples", numSamples);
			_feedrate.release_velo.normMB.temperature = T;
			_feedrate.release_velo.normMB.drift = D;
			_feedrate.release_velo.normMB.a_neg = a_neg;
			_feedrate.release_velo.normMB.a_pos = a_pos;
			_feedrate.release_velo.normMB.v_neg = v_neg;
			_feedrate.release_velo.normMB.v_pos = v_pos;
			_feedrate.release_velo.normMB.e_neg = e_neg;
			_feedrate.release_velo.normMB.e_pos = e_pos;
			_feedrate.release_velo.normMB.fpath = fpath;
			_feedrate.release_velo.numSamples = numSamples;

			// init random ins vector for released particles
			this->updateRandVecTrappedIns();

			// init random vectors generated for release velocities
			std::vector<double>& vx = _feedrate.release_velo.vx;
			std::vector<double>& vy = _feedrate.release_velo.vy;
			std::vector<double>& vz = _feedrate.release_velo.vz;
			update_velocity_vectors(_rnd.get(), numSamples, T, D, v_neg, e_neg, vx, vy, vz);
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

	// range|cid that is free of manipulation from MettDeamon
	{
		double ymin, ymax;
		ymin = ymax = 0.;
		uint32_t cid_ub = 0;
		xmlconfig.getNodeValue("control/manipfree/ymin", ymin);
		xmlconfig.getNodeValue("control/manipfree/ymax", ymax);
		xmlconfig.getNodeValue("control/manipfree/cid_ub", cid_ub);
		_manipfree.ymin = ymin;
		_manipfree.ymax = ymax;
		_manipfree.cid_ub = cid_ub;
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

	// change identity of fixed (frozen) molecules by component ID
	if(xmlconfig.changecurrentnode("changes")) {
		uint8_t numChanges = 0;
		XMLfile::Query query = xmlconfig.query("change");
		numChanges = query.card();
		global_log->info() << "[MettDeamon] Number of fixed molecules components: " << (uint32_t)numChanges << endl;
		if(numChanges < 1) {
			global_log->error() << "[MettDeamon] No component change defined in XML-config file. Program exit ..." << endl;
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
		global_log->error() << "[MettDeamon] No component changes defined in XML-config file. Program exit ..." << endl;
		Simulation::exit(-1);
	}
}

void MettDeamon::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain)
{
	double cutoffRadius = global_simulation->getcutoffRadius();
	this->prepare_start(domainDecomp, particleContainer, cutoffRadius);
}
void MettDeamon::beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep)
{
	this->init_positionMap(particleContainer);
}
void MettDeamon::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep)
{
	double cutoffRadius = global_simulation->getcutoffRadius();
	this->preForce_action(particleContainer, cutoffRadius);
}
void MettDeamon::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep)
{
	this->postForce_action(particleContainer, domainDecomp);
}

void MettDeamon::findMaxMoleculeID(DomainDecompBase* domainDecomp)
{
	ParticleContainer* particleContainer = global_simulation->getMoleculeContainer();
	_nMaxMoleculeID.local = 0;

	// max molecule id in particle container (local)
	for(auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
	{
		uint64_t id = pit->getID();
		if(id > _nMaxMoleculeID.local)
			_nMaxMoleculeID.local = id;
	}

	// global max IDs
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(_nMaxMoleculeID.local);
	domainDecomp->collCommAllreduceSum();
	_nMaxMoleculeID.global = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
}

uint64_t MettDeamon::getnNumMoleculesDeleted2( DomainDecompBase* domainDecomposition)
{
	domainDecomposition->collCommInit(1);
	domainDecomposition->collCommAppendUnsLong(_nNumMoleculesTooFast.local);
	domainDecomposition->collCommAllreduceSum();
	_nNumMoleculesTooFast.global = domainDecomposition->collCommGetUnsLong();
	domainDecomposition->collCommFinalize();
	return _nNumMoleculesTooFast.global;
}

void MettDeamon::prepare_start(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer, double cutoffRadius)
{
	_feedrate.feed.actual = _feedrate.feed.init;
	_reservoir->readParticleData(domainDecomp, particleContainer);
	_dInvDensityArea = 1. / (_dAreaXZ * _reservoir->getDensity(0) );
	if(_reservoir->getDensity(0) < 1e-9) {
		global_log->warning() << "[MettDeamon] ERROR: Reservoir density too low, _reservoir->getDensity(0)="
							<< _reservoir->getDensity(0) << endl;
	}
	// Activate reservoir bin with respect to restart information
	if(_bIsRestart)
		this->initRestart();

	this->InitTransitionPlane(global_simulation->getDomain() );

	// find max molecule ID in particle container
	this->findMaxMoleculeID(domainDecomp);

	particleContainer->deleteOuterParticles();

	if(_bIsRestart)
		return;

	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();

	for(auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
	{
		double dY = pit->r(1);
		if(dY > _manipfree.ymin && dY < _manipfree.ymax)
			continue;

		uint32_t cid_ub = pit->componentid()+1;
		if(cid_ub == _manipfree.cid_ub)
			continue;

		bool IsBehindTransitionPlane = this->IsBehindTransitionPlane(dY);
		if(not IsBehindTransitionPlane) {
			uint32_t cid = pit->componentid();
			if(cid != _vecChangeCompIDsFreeze.at(cid)) {
				Component* compNew = &(ptrComps->at(_vecChangeCompIDsFreeze.at(cid) ) );
				pit->setComponent(compNew);
			}
		}
	}
}

void MettDeamon::init_positionMap(ParticleContainer* particleContainer)
{
	for(auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
	{
		uint64_t mid = pit->getID();
		uint32_t cid = pit->componentid();

		bool bIsTrappedMolecule = this->IsTrappedMolecule(cid);
		if(bIsTrappedMolecule)
		{
			//savevelo
			std::array<double,10> pos;
			pos.at(0) = pit->r(0);
			pos.at(1) = pit->r(1);
			pos.at(2) = pit->r(2);
			pos.at(3) = pit->v(0);
			pos.at(4) = pit->v(1);
			pos.at(5) = pit->v(2);
			Quaternion q = pit->q();
			pos.at(6) = q.qw();
			pos.at(7) = q.qx();
			pos.at(8) = q.qy();
			pos.at(9) = q.qz();
			_storePosition[pit->getID()] = pos;
		}
	}
}

bool MettDeamon::IsInsideOuterReservoirSlab(const double& dPosY, const double& dBoxY)
{
	bool bRet = false;
	if(MD_LEFT_TO_RIGHT == _nMovingDirection)
		bRet = (dPosY < _reservoir->getBinWidth() );
	else if(MD_RIGHT_TO_LEFT == _nMovingDirection)
		bRet = (dPosY > (dBoxY - _reservoir->getBinWidth() ) );
	return bRet;
}

void MettDeamon::releaseTrappedMolecule(Molecule* mol, bool& bDeleteParticle)
{
	bDeleteParticle = false;
	uint16_t cid_zb = mol->componentid();
	if(not this->IsTrappedMolecule(cid_zb) || not this->IsBehindTransitionPlane(mol->r(1) ))
		return;

	// delete element from map to save memory
	size_t numDeleted = _storePosition.erase(mol->getID() );
	// only release particles with respect to parameter a_neg from normMB
	// they keep their component (trapped particle) although passing transition plane
	// feature DensityControl must delete them directly
	if(MD_RIGHT_TO_LEFT == _nMovingDirection) {
		// refill rand insertion vector to only release part of trapped particles, if no values left
		if(_feedrate.vec_rand_ins.empty())
			this->updateRandVecTrappedIns();

		int nDoInsert = _feedrate.vec_rand_ins.back();
		_feedrate.vec_rand_ins.pop_back();
		if(0 == nDoInsert)
		{
			mol->setv(1, 0.5);
			bDeleteParticle = true;
			_released.deleted.local++;
			return;
		}
	}

	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
	Component* compNew = &(ptrComps->at(_vecChangeCompIDsUnfreeze.at(cid_zb) ) );
	mol->setComponent(compNew);

	this->resetVelocity(mol);
	if(RVM_FIX_VALUE == _feedrate.release_velo.method)
		mol->setv(1, _feedrate.release_velo.fix_value);
	else if(RVM_ADD_FIX_VALUE == _feedrate.release_velo.method)
		mol->setv(1, mol->v(1) + _feedrate.release_velo.fix_value);
	else if(RVM_NORM_DISTR == _feedrate.release_velo.method) {
		float vxz = _norm.vxz.front();
		mol->setv(0, vxz);
		_norm.vxz.push_back(vxz);
		_norm.vxz.pop_front();
		vxz = _norm.vxz.front();
		mol->setv(2, vxz);
		_norm.vxz.push_back(vxz);
		_norm.vxz.pop_front();
		float vy = _norm.vy.front();
		mol->setv(1, vy);
		_norm.vy.push_back(vy);
		_norm.vy.pop_front();
	}
	else if(RVM_NORM_DISTR_GENERATOR == _feedrate.release_velo.method) {
		double T, D;  // T:temperature, D:Drift
		uint64_t N;  // numSamples
		T = _feedrate.release_velo.normMB.temperature;
		D = _feedrate.release_velo.normMB.drift;
		N = _feedrate.release_velo.numSamples;
		const double& v_neg = _feedrate.release_velo.normMB.v_neg;
		const double& e_neg = _feedrate.release_velo.normMB.e_neg;
		std::vector<double>& vx = _feedrate.release_velo.vx;
		std::vector<double>& vy = _feedrate.release_velo.vy;
		std::vector<double>& vz = _feedrate.release_velo.vz;

		if(vx.empty())
			update_velocity_vectors(_rnd.get(), N, T, D, v_neg, e_neg, vx, vy, vz);

		std::array<double,3> v;
		v[0] = vx.back(); vx.pop_back();
		v[1] = vy.back(); vy.pop_back();
		v[2] = vz.back(); vz.pop_back();

		_released.log_v.push_back(v);
		for(uint16_t dim=0; dim<3; dim++)
			mol->setv(dim, v[dim]);
	}
	// count released molecules
	_released.count.local++;
}

void MettDeamon::resetPositionAndOrientation(Molecule* mol, const double& dBoxY)
{
	uint16_t cid_zb = mol->componentid();
	if(not this->IsTrappedMolecule(cid_zb) )
		return;

	std::map<unsigned long, std::array<double,10> >::iterator it;
	it = _storePosition.find(mol->getID() );
	if(it == _storePosition.end() )
		return;

	// x,z position: always reset
	mol->setr(0, it->second.at(0) );
	mol->setr(2, it->second.at(2) );

	// reset y-position
	if(MD_LEFT_TO_RIGHT == _nMovingDirection)
		mol->setr(1, it->second.at(1) + _feedrate.feed.actual);
	else if(MD_RIGHT_TO_LEFT == _nMovingDirection)
		mol->setr(1, it->second.at(1) - _feedrate.feed.actual);

	// reset quaternion (orientation)
	Quaternion q(it->second.at(6),
			it->second.at(7),
			it->second.at(8),
			it->second.at(9) );
	mol->setq(q);
}

void MettDeamon::resetVelocity(Molecule* mol)
{
	uint16_t cid_zb = mol->componentid();
	if(not this->IsTrappedMolecule(cid_zb) )
		return;

	std::map<unsigned long, std::array<double,10> >::iterator it;
	it = _storePosition.find(mol->getID() );
	if(it == _storePosition.end() )
		return;

	// reset vx,vy,vz
	mol->setv(0, it->second.at(3) );
	mol->setv(1, it->second.at(4) );
	mol->setv(2, it->second.at(5) );
}

void MettDeamon::preForce_action(ParticleContainer* particleContainer, double cutoffRadius)
{
	double dBoxY = global_simulation->getDomain()->getGlobalLength(1);
	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();

	particleContainer->updateMoleculeCaches();

	for(auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
	{
		uint32_t cid = pit->componentid();
		uint32_t cid_ub = cid+1;
		double dY = pit->r(1);
		if(dY > _manipfree.ymin && dY < _manipfree.ymax)
			continue;

		if(cid_ub == _manipfree.cid_ub)
			continue;

		bool bIsTrappedMolecule = this->IsTrappedMolecule(cid);
		bool IsBehindTransitionPlane = this->IsBehindTransitionPlane(dY);

		// release trapped molecule
		bool bDeleteParticle = false;
		this->releaseTrappedMolecule( &(*pit), bDeleteParticle );

		if(bDeleteParticle)
		{
			particleContainer->deleteMolecule(pit, false);
			continue;
		}

		// reset position and orientation of fixed molecules
		this->resetPositionAndOrientation( &(*pit), dBoxY);

		if(bIsTrappedMolecule) {
			pit->setD(0, 0.);
			pit->setD(1, 0.);
			pit->setD(2, 0.);
		}
		this->resetVelocity( &(*pit) );

	}  // loop over molecules

	// Nothing more to do in case of empty reservoir
	if(RRM_EMPTY == _reservoir->getReadMethod() )
		return;

	_feedrate.feed.sum += _feedrate.feed.actual;
	if (_feedrate.feed.sum >= _reservoir->getBinWidth())
	{
		// insert actual reservoir slab / activate next reservoir slab
		DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
		this->findMaxMoleculeID(&domainDecomp);
		this->InsertReservoirSlab(particleContainer);
	}

	// log count of released molecules
	this->logReleased();

	// write released velocities
	this->logReleasedVelocities();
}
void MettDeamon::postForce_action(ParticleContainer* particleContainer, DomainDecompBase* domainDecomposition)
{
	// Nothing to do in case of empty reservoir
	if(RRM_EMPTY == _reservoir->getReadMethod() )
		return;

	unsigned long nNumMoleculesLocal = 0;
	unsigned long nNumMoleculesGlobal = 0;

	for(auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
	{
		uint32_t cid_zb = pit->componentid();
		uint32_t cid_ub = cid_zb+1;
		double dY = pit->r(1);
		if(dY > _manipfree.ymin && dY < _manipfree.ymax)
			continue;

		if(cid_ub == _manipfree.cid_ub)
			continue;

		bool bIsTrappedMolecule = this->IsTrappedMolecule(cid_zb);

		if(bIsTrappedMolecule) {
			pit->setD(0, 0.);
			pit->setD(1, 0.);
			pit->setD(2, 0.);

			this->resetVelocity( &(*pit) );
		}

	}  // loop over molecules

	nNumMoleculesLocal = particleContainer->getNumberOfParticles();

	// Update feedrate
	if( (FRM_DIRECTED == _nFeedRateMethod) && (global_simulation->getSimulationStep() % _nUpdateFreq == 0) )
	{
		// update global number of particles / calc global number of deleted particles
		domainDecomposition->collCommInit(5);
		domainDecomposition->collCommAppendUnsLong(nNumMoleculesLocal);
		domainDecomposition->collCommAppendUnsLong(_feedrate.numMolecules.inserted.local);
		domainDecomposition->collCommAppendUnsLong(_feedrate.numMolecules.deleted.local);
		domainDecomposition->collCommAppendUnsLong(_feedrate.numMolecules.changed_to.local);
		domainDecomposition->collCommAppendUnsLong(_feedrate.numMolecules.changed_from.local);
		domainDecomposition->collCommAllreduceSum();
		nNumMoleculesGlobal = domainDecomposition->collCommGetUnsLong();
		_feedrate.numMolecules.inserted.global = domainDecomposition->collCommGetUnsLong();
		_feedrate.numMolecules.deleted.global = domainDecomposition->collCommGetUnsLong();
		_feedrate.numMolecules.changed_to.global = domainDecomposition->collCommGetUnsLong();
		_feedrate.numMolecules.changed_from.global = domainDecomposition->collCommGetUnsLong();
		domainDecomposition->collCommFinalize();
		_nNumMoleculesDeletedGlobalAlltime += _feedrate.numMolecules.deleted.global;
		_feedrate.numMolecules.inserted.local = 0;
		_feedrate.numMolecules.deleted.local = 0;
		_feedrate.numMolecules.changed_to.local = 0;
		_feedrate.numMolecules.changed_from.local = 0;

		// update sum and summation list
		int64_t numMolsDeletedOrChanged = 0;
		if(FRM_DELETED_MOLECULES == _nFeedRateMethod)
		{
			numMolsDeletedOrChanged += _feedrate.numMolecules.deleted.global;
			numMolsDeletedOrChanged -= _feedrate.numMolecules.inserted.global;
			numMolsDeletedOrChanged += _feedrate.numMolecules.changed_from.global;
			numMolsDeletedOrChanged -= _feedrate.numMolecules.changed_to.global;
		}
		else if(FRM_CHANGED_MOLECULES == _nFeedRateMethod)
		{
			numMolsDeletedOrChanged += _feedrate.numMolecules.changed_from.global;
			numMolsDeletedOrChanged -= _feedrate.numMolecules.changed_to.global;
		}
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
		global_log->debug() << "_numDeletedMolsSum = " << _numDeletedMolsSum << endl;
		global_log->debug() << "_dDeletedMolsPerTimestep = " << _dDeletedMolsPerTimestep << endl;
		global_log->debug() << "_feedrate.feed.actual = " << _feedrate.feed.actual << endl;
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

	// log feed rate
	this->logFeedrate();
}

void MettDeamon::writeRestartfile()
{
	uint64_t simstep = global_simulation->getSimulationStep();
	if(0 != simstep % _nWriteFreqRestart)
		return;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

	if(domainDecomp.getRank() != 0)
		return;

	// init restart file
	if(not _bInitRestartLog)
	{
		std::stringstream fnamestream;
		fnamestream << "MettDeamonRestart_movdir-" << (uint32_t)_nMovingDirection << ".dat";
		std::ofstream ofs(fnamestream.str().c_str(), std::ios::out);
		std::stringstream outputstream;
		outputstream << "     simstep" << "   slabIndex" << "                  deltaY" << std::endl;
		ofs << outputstream.str();
		ofs.close();
		_bInitRestartLog = true;
	}

	std::stringstream fnamestream;
	fnamestream << "MettDeamonRestart_movdir-" << (uint32_t)_nMovingDirection << ".dat";
	std::ofstream ofs(fnamestream.str().c_str(), std::ios::app);
	std::stringstream outputstream;

	outputstream << setw(12) << simstep << setw(12) << _reservoir->getActualBinIndex();
	outputstream << FORMAT_SCI_MAX_DIGITS << _feedrate.feed.sum << std::endl;

	ofs << outputstream.str();
	ofs.close();

	// write restart info in XML format
	{
		std::stringstream fnamestream;
		fnamestream << "MettDeamonRestart_movdir-" << (uint32_t)_nMovingDirection << "_TS" << fill_width('0', 9) << simstep << ".xml";
		std::ofstream ofs(fnamestream.str().c_str(), std::ios::out);
		std::stringstream outputstream;
		ofs << "<?xml version='1.0' encoding='UTF-8'?>" << endl;
		ofs << "<restart>" << endl;
		ofs << "\t<binindex>" << _reservoir->getActualBinIndex() << "</binindex>" << endl;
		ios::fmtflags f( ofs.flags() );
		ofs << "\t<deltaY>" << FORMAT_SCI_MAX_DIGITS_WIDTH_21 << _feedrate.feed.sum << "</deltaY>" << endl;
		ofs.flags(f);  // restore default format flags
		ofs << "</restart>" << endl;

		ofs << outputstream.str();
		ofs.close();
	}
}

void MettDeamon::logFeedrate()
{
	if(0 != global_simulation->getSimulationStep() % _feedrate.log_freq)
		return;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

	if(domainDecomp.getRank() != 0)
		return;

	// init feedrate log file
	if(not _bInitFeedrateLog)
	{
		std::stringstream fnamestream;
		fnamestream << "MettDeamon_feedrate_movdir-" << (uint32_t)_nMovingDirection << ".dat";
		std::ofstream ofs(fnamestream.str().c_str(), std::ios::out);
		std::stringstream outputstream;
		outputstream << "     simstep" << "                feedrate" << std::endl;
		ofs << outputstream.str();
		ofs.close();
		_bInitFeedrateLog = true;
	}

	std::stringstream fnamestream;
	fnamestream << "MettDeamon_feedrate_movdir-" << (uint32_t)_nMovingDirection << ".dat";
	std::ofstream ofs(fnamestream.str().c_str(), std::ios::app);
	std::stringstream outputstream;

	outputstream << setw(12) << global_simulation->getSimulationStep();
	outputstream << FORMAT_SCI_MAX_DIGITS << _feedrate.feed.actual << std::endl;

	ofs << outputstream.str();
	ofs.close();
}

void MettDeamon::logReleased()
{
	if(0 != global_simulation->getSimulationStep() % _feedrate.log_freq)
		return;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	domainDecomp.collCommInit(2);
	domainDecomp.collCommAppendUnsLong(_released.count.local);
	domainDecomp.collCommAppendUnsLong(_released.deleted.local);
	domainDecomp.collCommAllreduceSum();
	_released.count.global = domainDecomp.collCommGetUnsLong();
	_released.deleted.global = domainDecomp.collCommGetUnsLong();
	domainDecomp.collCommFinalize();

	// reset local values
	_released.count.local = 0;
	_released.deleted.local = 0;

	if(domainDecomp.getRank() != 0)
		return;

	// init released count log file
	if(not _released.init_file)
	{
		std::stringstream fnamestream;
		fnamestream << "MettDeamon_released_movdir-" << (uint32_t)_nMovingDirection << ".dat";
		std::ofstream ofs(fnamestream.str().c_str(), std::ios::out);
		std::stringstream outputstream;
		outputstream << "     simstep" << "       count" << "     deleted" << std::endl;
		ofs << outputstream.str();
		ofs.close();
		_released.init_file = true;
	}

	std::stringstream fnamestream;
	fnamestream << "MettDeamon_released_movdir-" << (uint32_t)_nMovingDirection << ".dat";
	std::ofstream ofs(fnamestream.str().c_str(), std::ios::app);
	std::stringstream outputstream;

	outputstream << setw(12) << global_simulation->getSimulationStep() << setw(12) << _released.count.global << setw(12) << _released.deleted.global << std::endl;

	ofs << outputstream.str();
	ofs.close();
}

void MettDeamon::logReleasedVelocities()
{
	if(_released.log_v.empty() )
		return;

	uint64_t simstep = global_simulation->getSimulationStep();
	if(0 != (simstep % _released.log_freq_vel) )
		return;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	int nRank = domainDecomp.getRank();

	// construct filename
	std::stringstream fnamestream;
	fnamestream << "MettDeamon_released_vel_movdir-" << (uint32_t)_nMovingDirection << "_TS" << fill_width('0', 9) << simstep << "_p" << nRank << ".dat";

	std::ofstream ofs(fnamestream.str().c_str(), std::ios::out);
	std::stringstream outputstream;
	outputstream << "                      vx" << "                      vy" << "                      vz" << std::endl;

	for(auto vi:_released.log_v)
	{
		outputstream << FORMAT_SCI_MAX_DIGITS << vi[0] << FORMAT_SCI_MAX_DIGITS << vi[1] << FORMAT_SCI_MAX_DIGITS << vi[2] << std::endl;
	}
	_released.log_v.clear();

	ofs << outputstream.str();
	ofs.close();
}

void MettDeamon::calcDeltaY()
{
	_feedrate.feed.actual = _dDeletedMolsPerTimestep * _dInvDensityArea;
	if(MD_LEFT_TO_RIGHT == _nMovingDirection && _feedrate.feed.actual < 0.)
		_feedrate.feed.actual = 0.;
	else if (MD_RIGHT_TO_LEFT == _nMovingDirection && _feedrate.feed.actual < 0.)  // > 0.
		_feedrate.feed.actual = 0.;
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
		_feedrate.feed.actual = 0.;
	else
		_feedrate.feed.actual = dDensityDelta/_reservoir->getDensity(0)*dInvNumVals*_dVolumeCV/_dAreaXZ;
}

void MettDeamon::InitTransitionPlane(Domain* domain)
{
	double dBoxY = domain->getGlobalLength(1);
	if(MD_LEFT_TO_RIGHT == _nMovingDirection)
		_dTransitionPlanePosY = 2*_reservoir->getBinWidth();
	else
		_dTransitionPlanePosY = dBoxY - 2*_reservoir->getBinWidth();
}

void MettDeamon::getAvailableParticleIDs(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		CommVar<std::vector<uint64_t> >& particleIDs_available, const CommVar<uint64_t>& numParticleIDs)
{
	int nRank = domainDecomp->getRank();
	int numProcs = domainDecomp->getNumProcs();
	Domain* domain = global_simulation->getDomain();
	CommVar<std::vector<uint32_t> > particleIDs_assigned;
	CommVar<uint64_t> maxID;
	CommVar<uint64_t> numMolecules;
	domain->updateMaxMoleculeID(particleContainer, domainDecomp);
	maxID = domain->getMaxMoleculeID();
	domain->updateglobalNumMolecules(particleContainer, domainDecomp);
	numMolecules.global = domain->getglobalNumMolecules();
	global_log->debug() << "[" << nRank << "]: maxID.local, maxID.global=" << maxID.local << ", " << maxID.global << endl;
	uint64_t numMoleculesAfterInsertion = numMolecules.global + numParticleIDs.global;
	uint64_t numIDs;
	if(maxID.global >= numMoleculesAfterInsertion)
		numIDs = maxID.global + 1;
	else
		numIDs = numMoleculesAfterInsertion + 1;
	global_log->debug() << "[" << nRank << "]: numMoleculesAfterInsertion, numMolecules.global=" << numMoleculesAfterInsertion << ", " << numMolecules.global << endl;

	std::vector<uint32_t>& vl = particleIDs_assigned.local;
	vl.resize(numIDs); std::fill(vl.begin(), vl.end(), 0);
	std::vector<uint32_t>& vg = particleIDs_assigned.global;
	vg.resize(numIDs); std::fill(vg.begin(), vg.end(), 0);

	// mark 0:available | 1:assigned particle IDs, 0:
	for(auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
		uint64_t pid = pit->getID();
		vl.at(pid) = 1;
	}
#ifdef ENABLE_MPI
	MPI_Allreduce(vl.data(), vg.data(), vg.size(), MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
#else
	for(auto ii=0; ii<vg.size(); ++ii)
		vg.at(ii) = vl.at(ii);
#endif

	global_log->debug() << "[" << nRank << "]: assigned.local, assigned.global=" << particleIDs_assigned.local.size() << ", " << particleIDs_assigned.global.size() << endl;

	// check for available particle IDs
	for(uint64_t pid=1; pid<vg.size() && particleIDs_available.global.size()<numParticleIDs.global; ++pid) {
		if(0 == vg.at(pid))
			particleIDs_available.global.push_back(pid);
	}
	particleIDs_available.local.resize(numParticleIDs.local);
	global_log->debug() << "[" << nRank << "]: avail.local, avail.global=" << particleIDs_available.local.size() << ", " << particleIDs_available.global.size() << endl;

#ifdef ENABLE_MPI
	// gather displacement (displs)
	std::vector<int32_t> counts, displs;
	counts.resize(numProcs); displs.resize(numProcs);
	{
	int32_t sendbuf = (int32_t)(particleIDs_available.local.size() );
	int32_t* recvbuf = NULL;
	if(0 == nRank)
		recvbuf = counts.data();
	MPI_Gather(&sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	// scatter available particle IDs
	{
	displs.at(0) = 0;
	for (auto ri=1; ri<numProcs; ++ri)
		displs.at(ri) = displs.at(ri-1) + counts.at(ri-1);

	uint64_t* sendbuf = NULL;
	if(0 == nRank)
		sendbuf = particleIDs_available.global.data();
	MPI_Scatterv(sendbuf,
			counts.data(),
			displs.data(),
			MPI_UNSIGNED_LONG,
			particleIDs_available.local.data(),
			particleIDs_available.local.size(),
			MPI_UNSIGNED_LONG,
			0,
			MPI_COMM_WORLD);
	}
#else
	particleIDs_available.local = particleIDs_available.global;
#endif
}

void MettDeamon::updateReservoir(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer)
{
	_reservoir->updateParticleData(domainDecomp, particleContainer);
}

void MettDeamon::InsertReservoirSlab(ParticleContainer* particleContainer)
{
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	this->updateReservoir(&domainDecomp, particleContainer);
	
	int nRank = domainDecomp.getRank();
	int numProcs = domainDecomp.getNumProcs();
	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
	std::vector<Molecule>& currentReservoirSlab = _reservoir->getParticlesActualBin();
	_reservoir->printBinQueueInfo();

	CommVar<uint64_t> numParticlesCurrentSlab;
	numParticlesCurrentSlab.local = currentReservoirSlab.size();
	// calc global values
	domainDecomp.collCommInit(1);
	domainDecomp.collCommAppendUnsLong(numParticlesCurrentSlab.local);
	domainDecomp.collCommAllreduceSum();
	numParticlesCurrentSlab.global = domainDecomp.collCommGetUnsLong();
	domainDecomp.collCommFinalize();

	// get available particle IDs
	CommVar<std::vector<uint64_t> > particleIDs_available;
	this->getAvailableParticleIDs(particleContainer, &domainDecomp, particleIDs_available, numParticlesCurrentSlab);

	// reduce reservoir density to percentage
	CommVar<uint64_t> numAdded;
	numAdded.local = 0;
	double percent = _reservoir->GetInsPercent();
	std::vector<int> v;
	create_rand_vec_ones(numParticlesCurrentSlab.local, percent, v);
	int64_t index = -1;

	for(auto mi : currentReservoirSlab)
	{
		index++;
		uint64_t id = mi.getID();
		uint32_t cid = mi.componentid();
		Component* compNew;
		if(_nMovingDirection==1 || _nMovingDirection==2)
			compNew = &(ptrComps->at(_vecChangeCompIDsFreeze.at(cid) ) );
		else
			compNew = &(ptrComps->at(3));

		mi.setid(particleIDs_available.local.at(index) );
		mi.setComponent(compNew);
		mi.setr(1, mi.r(1) + _feedrate.feed.sum - _reservoir->getBinWidth() );
		particleContainer->addParticle(mi);
		numAdded.local++;
	}
	_feedrate.feed.sum -= _reservoir->getBinWidth();  // reset feed sum
	if(not _reservoir->nextBin(_nMaxMoleculeID.global) ) {
		global_log->error() << "[MettDeamon] Failed to activate new bin of particle Reservoir's BinQueue => Program exit." << endl;
		Simulation::exit(-1);
	}
	global_log->debug() << "[" << nRank << "]: ADDED " << numAdded.local << "/" << numParticlesCurrentSlab.local << " particles (" << numAdded.local/static_cast<float>(numParticlesCurrentSlab.local)*100 << ")%." << endl;
	// calc global values
	domainDecomp.collCommInit(1);
	domainDecomp.collCommAppendUnsLong(numAdded.local);
	domainDecomp.collCommAllreduceSum();
	numAdded.global = domainDecomp.collCommGetUnsLong();
	domainDecomp.collCommFinalize();

	if(0 == nRank)
		global_log->debug() << "[" << nRank << "]: ADDED " << numAdded.global << "/" << numParticlesCurrentSlab.global << " particles (" << numAdded.global/static_cast<float>(numParticlesCurrentSlab.global)*100 << ")%." << endl;
}

void MettDeamon::initRestart()
{
	bool bRet = _reservoir->activateBin(_restartInfo.nBindindex);
	if(not bRet)
	{
		global_log->info() << "[MettDeamon] Failed to activate reservoir bin after restart! Program exit ... " << endl;
		Simulation::exit(-1);
	}
	_feedrate.feed.sum = _restartInfo.dYsum;
}

void MettDeamon::readNormDistr()
{
	struct {
		std::ifstream vxz;
		std::ifstream vy;
	} ifs{};
	ifs.vxz.open(_norm.fname.vxz, std::ios::in);
	ifs.vy.open(_norm.fname.vy, std::ios::in);

	//check to see that the file was opened correctly:
	if (!ifs.vxz.is_open() || !ifs.vy.is_open() ) {
		std::cerr << "[MettDeamon] There was a problem opening the input file!\n";
		Simulation::exit(-1);//exit or do additional error checking
	}

	double dVal = 0.0;
	//keep storing values from the text file so long as data exists:
	while (ifs.vxz >> dVal) {
		_norm.vxz.push_back(dVal);
	}
	while (ifs.vy >> dVal) {
		if(MD_LEFT_TO_RIGHT == _nMovingDirection)
			_norm.vy.push_back( abs(dVal) );
		else if (MD_RIGHT_TO_LEFT == _nMovingDirection)
			_norm.vy.push_back( abs(dVal) * (-1.) );
		else
			Simulation::exit(-1);
	}
	// close files
	ifs.vxz.close();
	ifs.vy.close();
}

void MettDeamon::updateRandVecTrappedIns()
{
	create_rand_vec_ones(100, _feedrate.release_velo.normMB.a_neg, _feedrate.vec_rand_ins);
	global_log->debug() << "_feedrate.vec_rand_ins: ";
	int nSum = 0;
	for(auto vi:_feedrate.vec_rand_ins) {
		global_log->debug() << vi << ",";
		nSum += vi;
	}
	global_log->debug() << "sum=" << nSum << endl;
}

// class Reservoir
Reservoir::Reservoir(MettDeamon* parent) :
	_parent(parent),
	_moleculeDataReader(nullptr),
	_binQueue(nullptr),
	_numMoleculesRead(0),
	_nMaxMoleculeID(0),
	_nMoleculeFormat(ICRVQD),
	_nReadMethod(RRM_UNKNOWN),
	_dReadWidthY(0.0),
	_dBinWidthInit(0.0),
	_dBinWidth(0.0)
{
	// init filepath struct
	_filepath.data = _filepath.header = "unknown";

	// allocate BinQueue
	_binQueue.reset(new BinQueue());

	// init identity change vector
	uint16_t nNumComponents = global_simulation->getEnsemble()->getComponents()->size();
	_vecChangeCompIDs.resize(nNumComponents);
	std::iota (std::begin(_vecChangeCompIDs), std::end(_vecChangeCompIDs), 0);

	// init density vector
	_density.resize(nNumComponents+1);  // 0: total density
}

Reservoir::~Reservoir() = default;

void Reservoir::readXML(XMLfileUnits& xmlconfig)
{
	// update BinQueue before inserting new Reservoir slab
	_bUpdateBinQueue = false;
	xmlconfig.getNodeValue("@update", _bUpdateBinQueue);
	
	std::string strType = "unknown";
	bool bRet1 = xmlconfig.getNodeValue("file@type", strType);
	bool bRet2 = xmlconfig.getNodeValue("width", _dReadWidthY);
	xmlconfig.getNodeValue("binwidth", _dBinWidthInit);
	_dInsPercent = 1.0;
	xmlconfig.getNodeValue("ins_percent", _dInsPercent);
	_nReadMethod = RRM_UNKNOWN;
	if(bRet1)
	{
		if("ASCII" == strType) {
			_nReadMethod = RRM_READ_FROM_FILE;
			xmlconfig.getNodeValue("file", _filepath.data);
			_filepath.header = _filepath.data;
		}
		else if("binary" == strType) {
			_nReadMethod = RRM_READ_FROM_FILE_BINARY;
			xmlconfig.getNodeValue("file/header", _filepath.header);
			xmlconfig.getNodeValue("file/data", _filepath.data);
		}
		else if("empty" == strType) {
			_nReadMethod = RRM_EMPTY;
		}
		else {
			global_log->error() << "[MettDeamon] Wrong file type='" << strType << "' specified. Programm exit ..." << endl;
			Simulation::exit(-1);
		}
	}
	else if(not bRet1 and (RRM_EMPTY != _nReadMethod) )
		_nReadMethod = RRM_READ_FROM_MEMORY;
	else
		_nReadMethod = RRM_AMBIGUOUS;

	// Possibly change component IDs
	if(xmlconfig.changecurrentnode("changes")) {
		uint8_t numChanges = 0;
		XMLfile::Query query = xmlconfig.query("change");
		numChanges = query.card();
		if(numChanges < 1) {
			global_log->error() << "[MettDeamon] No component change defined in XML-config file. Program exit ..." << endl;
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

void Reservoir::readParticleData(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer)
{
	switch(_nReadMethod)
	{
	case RRM_READ_FROM_MEMORY:
	case RRM_EMPTY:
		this->readFromMemory(domainDecomp, particleContainer);
		break;
	case RRM_READ_FROM_FILE:
		this->readFromFile(domainDecomp, particleContainer);
		break;
	case RRM_READ_FROM_FILE_BINARY:
		this->readFromFileBinary(domainDecomp, particleContainer);
		break;
	case RRM_UNKNOWN:
	case RRM_AMBIGUOUS:
	default:
		global_log->error() << "[MettDeamon] Unknown (or ambiguous) method to read reservoir for feature MettDeamon. Program exit ..." << endl;
		Simulation::exit(-1);
	}

	// sort particles into bins
	this->sortParticlesToBins(domainDecomp, particleContainer);

	// volume, densities
	this->calcPartialDensities(domainDecomp);
}

void Reservoir::updateParticleData(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer)
{
	if (not _bUpdateBinQueue) {
		return;
	}
	Domain* domain = global_simulation->getDomain();
#ifndef ENABLE_MPI
	return;
#else

#define PARTICLE_BUFFER_SIZE  (16*1024)
	int ownRank = domainDecomp->getRank();

	/* distribute molecules to other MPI processes */
	for (int rank = 0; rank < domainDecomp->getNumProcs(); rank++) {
		unsigned long num_particles = _particleVector.size();
		MPI_CHECK( MPI_Bcast(&num_particles, 1, MPI_UNSIGNED_LONG, rank, domainDecomp->getCommunicator()) );

		ParticleData particle_buff[PARTICLE_BUFFER_SIZE];
		int particle_buff_pos = 0;
		MPI_Datatype mpi_Particle;
		ParticleData::getMPIType(mpi_Particle);

		if (rank == ownRank) {
			for(unsigned long i = 0; i < num_particles; ++i) {
				ParticleData::MoleculeToParticleData(particle_buff[particle_buff_pos], _particleVector[i]);
				particle_buff_pos++;
				if ((particle_buff_pos >= PARTICLE_BUFFER_SIZE) || (i == num_particles - 1)) {
					global_log->debug() << "broadcasting(sending) particles" << endl;
					MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, rank, domainDecomp->getCommunicator());
					particle_buff_pos = 0;
				}
			}
		} else {
			uint64_t numParticlesAdd = 0;
			for(unsigned long i = 0; i < num_particles; ++i) {
				if(i % PARTICLE_BUFFER_SIZE == 0) {
					global_log->debug() << "broadcasting(receiving) particles" << endl;
					MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, rank, domainDecomp->getCommunicator());
					particle_buff_pos = 0;
				}
				Molecule mol;
				ParticleData::ParticleDataToMolecule(particle_buff[particle_buff_pos], mol);
				particle_buff_pos++;
				
				bool bIsRelevant = this->isRelevant(domainDecomp, domain, mol);
				if (bIsRelevant) {
					_particleVector.push_back(mol);
					numParticlesAdd++;
				}
			}
			if(numParticlesAdd > 0) {
				global_log->debug() << "Rank " << ownRank << " received " << numParticlesAdd << " particles from rank " << rank << "." << endl;
			}
		}
		global_log->debug() << "broadcasting(sending/receiving) particles complete" << endl;
	}

	// delete particles out of bounding box
	uint64_t numParticlesOld = _particleVector.size();
	std::vector<Molecule> particleVectorTmp;
	for (auto mol : _particleVector) {
		bool bIsRelevant = this->isRelevant(domainDecomp, domain, mol);
		if (bIsRelevant) {
			particleVectorTmp.push_back(mol);
		}
	}
	_particleVector.resize(particleVectorTmp.size());
	_particleVector = particleVectorTmp;
	if(_particleVector.size() < numParticlesOld) {
		global_log->debug() << "Rank " << ownRank << " deleted " << numParticlesOld - _particleVector.size() << " particles from particle vector." << endl;
	}
	// Refresh BinQueue
	uint32_t actual = _binQueue->getActualBinIndex();
	this->clearBinQueue();
	this->sortParticlesToBins(domainDecomp, particleContainer);
	_binQueue->activateBin(actual);
#endif  //ENABLE_MPI
}

void Reservoir::sortParticlesToBins(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer)
{
	Domain* domain = global_simulation->getDomain();

	uint32_t numBins = _box.length.at(1) / _dBinWidthInit;
	_dBinWidth = _box.length.at(1) / (double)(numBins);
	global_log->debug() << "_arrBoxLength[1]="<<_box.length.at(1)<<endl;
	global_log->debug() << "_dBinWidthInit="<<_dBinWidthInit<<endl;
	global_log->debug() << "_numBins="<<numBins<<endl;
	global_log->debug() << "_particleVector.size()=" << _particleVector.size() << endl;

	global_log->debug() << "bbMin ="
			<< domainDecomp->getBoundingBoxMin(0, domain) << ", "
			<< domainDecomp->getBoundingBoxMin(1, domain) << ", "
			<< domainDecomp->getBoundingBoxMin(2, domain) << "; bbMax = "
			<< domainDecomp->getBoundingBoxMax(0, domain) << ", "
			<< domainDecomp->getBoundingBoxMax(1, domain) << ", "
			<< domainDecomp->getBoundingBoxMax(2, domain) << endl;
	std::vector< std::vector<Molecule> > binVector;
	binVector.resize(numBins);
	uint32_t nBinIndex;
	for(auto mol:_particleVector)  // copy of the molecule is needed, as we modify it.
	{
		// possibly change component IDs
		this->changeComponentID(mol, mol.componentid() );
		double y = mol.r(1);
		nBinIndex = floor(y / _dBinWidth);
		global_log->debug() << "y="<<y<<", nBinIndex="<<nBinIndex<<", _binVector.size()="<<binVector.size()<<endl;
		mardyn_assert(nBinIndex < binVector.size() );
		switch(_parent->getMovingDirection() )
		{
		case MD_LEFT_TO_RIGHT:
			mol.setr(1, y - nBinIndex*_dBinWidth);  // positions in slabs related to origin (x,y,z) == (0,0,0)
			break;
		case MD_RIGHT_TO_LEFT:
			mol.setr(1, y - nBinIndex*_dBinWidth + (domain->getGlobalLength(1) - _dBinWidth) );  // positions in slabs related to origin (x,y,z) == (0,0,0)
			break;
		}
		// check if molecule is in bounding box of the process domain
		bool bIsInsideBB = domainDecomp->procOwnsPos(mol.r(0), mol.r(1), mol.r(2), domain);
		bool bIsInsidePC = particleContainer->isInBoundingBox(mol.r_arr().data());
		if(bIsInsideBB != bIsInsidePC)
			global_log->debug() << "bIsInsideBB=" << bIsInsideBB << ", bIsInsidePC=" << bIsInsidePC << endl;
		if (bIsInsideBB)
			binVector.at(nBinIndex).push_back(mol);
	}

	// add bin particle vectors to bin queue
	switch(_parent->getMovingDirection() )
	{
		case MD_LEFT_TO_RIGHT:
			for (auto bit = binVector.rbegin(); bit != binVector.rend(); ++bit)
			{
				global_log->debug() << "(*bit).size()=" << (*bit).size() << endl;
				_binQueue->enque(*bit);
			}
			break;

		case MD_RIGHT_TO_LEFT:
			for(const auto& bin:binVector)
			{
				global_log->debug() << "bin.size()=" << bin.size() << endl;
				_binQueue->enque(bin);
			}
			break;
	}
}

void Reservoir::readFromMemory(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer)
{
	Domain* domain = global_simulation->getDomain();

	_box.length.at(0) = domain->getGlobalLength(0);
	_box.length.at(1) = _dReadWidthY;
	_box.length.at(2) = domain->getGlobalLength(2);

	if(RRM_EMPTY == _nReadMethod)
		return;

	for(auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
	{
		Molecule mol(*pit);
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

void Reservoir::readFromFile(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer)
{
	Domain* domain = global_simulation->getDomain();
	std::ifstream ifs;
	global_log->info() << "[MettDeamon] Opening Reservoirfile " << _filepath.data << endl;
	ifs.open( _filepath.data.c_str() );
	if (!ifs.is_open()) {
		global_log->error() << "[MettDeamon] Could not open Mettdeamon Reservoirfile " << _filepath.data << endl;
		Simulation::exit(1);
	}
	global_log->info() << "[MettDeamon] Reading Mettdeamon Reservoirfile " << _filepath.data << endl;

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
			_box.length.at(0) = domain->getGlobalLength(0);
			_box.length.at(1) = Ylength;
			_box.length.at(2) = domain->getGlobalLength(2);
		}
	}

	if((token != "NumberOfMolecules") && (token != "N")) {
		global_log->error() << "[MettDeamon] Expected the token 'NumberOfMolecules (N)' instead of '" << token << "'" << endl;
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
			global_log->error() << "[MettDeamon] Unknown molecule format '" << ntypestring << "'" << endl;
			Simulation::exit(1);
		}
	} else {
		ifs.seekg(spos);
	}
	global_log->info() << " molecule format: " << ntypestring << endl;

	if( numcomponents < 1 ) {
		global_log->warning() << "[MettDeamon] No components defined! Setting up single one-centered LJ" << endl;
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
			global_log->error() << "[MettDeamon] Molecule id " << id << " has wrong componentid: " << componentid << ">" << numcomponents << endl;
			Simulation::exit(1);
		}
		componentid --; // TODO: Component IDs start with 0 in the program.
		Molecule mol = Molecule(i+1,&dcomponents[componentid],x,y,z,vx,vy,vz,q0,q1,q2,q3,Dx,Dy,Dz);

		bool bIsRelevant = this->isRelevant(domainDecomp, domain, mol);
		if (bIsRelevant) {
			_particleVector.push_back(mol);
		}

		// Print status message
		unsigned long iph = _numMoleculesRead / 100;
		if( iph != 0 && (i % iph) == 0 )
			global_log->info() << "[MettDeamon] Finished reading molecules: " << i/iph << "%\r" << flush;
	}

	ifs.close();
}

void Reservoir::readFromFileBinaryHeader()
{
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	XMLfileUnits inp(_filepath.header);

	if(not inp.changecurrentnode("/mardyn")) {
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
	this->setDensity(0, nNumMols / dVolume);

	if(not bInputOk)
	{
		global_log->error() << "[MettDeamon] Content of file: '" << _filepath.header << "' corrupted! Program exit ..." << endl;
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
		global_log->error() << "[MettDeamon] Not a valid molecule format: " << strMoleculeFormat << ", program exit ..." << endl;
		Simulation::exit(1);
	}
}

void Reservoir::readFromFileBinary(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer)
{
	Domain* domain = global_simulation->getDomain();
	global_log->info() << "[MettDeamon] Reservoir::readFromFileBinary(...)" << endl;
	// read header
	this->readFromFileBinaryHeader();

#ifdef ENABLE_MPI
	if(domainDecomp->getRank() == 0) {
#endif
	global_log->info() << "[MettDeamon] Opening phase space file " << _filepath.data << endl;
	std::ifstream ifs;
	ifs.open(_filepath.data.c_str(), ios::binary | ios::in);
	if (!ifs.is_open()) {
		global_log->error() << "[MettDeamon] Could not open phaseSpaceFile " << _filepath.data << endl;
		Simulation::exit(1);
	}

	global_log->info() << "[MettDeamon] Reading phase space file " << _filepath.data << endl;

	vector<Component>& components = *(_simulation.getEnsemble()->getComponents());

	// Select appropriate reader
	switch (_nMoleculeFormat) {
		case ICRVQD: _moleculeDataReader.reset(new MoleculeDataReaderICRVQD()); break;
		case ICRV: _moleculeDataReader.reset(new MoleculeDataReaderICRV()); break;
		case IRV: _moleculeDataReader.reset(new MoleculeDataReaderIRV()); break;
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
		std::vector<Molecule> particleVectorTmp;
		for(unsigned long i = 0; i < num_particles; ++i) {
			ParticleData::MoleculeToParticleData(particle_buff[particle_buff_pos], _particleVector[i]);
			particle_buff_pos++;
			if ((particle_buff_pos >= PARTICLE_BUFFER_SIZE) || (i == num_particles - 1)) {
				global_log->debug() << "broadcasting(sending) particles" << endl;
				MPI_Bcast(particle_buff, PARTICLE_BUFFER_SIZE, mpi_Particle, 0, domainDecomp->getCommunicator());
				particle_buff_pos = 0;
			}
			Molecule& mol = _particleVector[i];
			bool bIsRelevant = this->isRelevant(domainDecomp, domain, mol);
			if (bIsRelevant) {
				particleVectorTmp.push_back(mol);
			}
		}
		_particleVector.resize(particleVectorTmp.size());
		_particleVector = particleVectorTmp;
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
			
			bool bIsRelevant = this->isRelevant(domainDecomp, domain, mol);
			if (bIsRelevant) {
				_particleVector.push_back(mol);
			}
		}
	}
	global_log->debug() << "broadcasting(sending/receiving) particles complete" << endl;
#endif
}

void Reservoir::calcPartialDensities(DomainDecompBase* domainDecomp)
{
	// calc box volume
	_box.volume = 1.; for(uint8_t dim=0; dim<3; dim++) _box.volume *= _box.length[dim];
	double dInvVolume = 1./_box.volume;

	// count particles of each component
	// TODO: not nice that in case of RRM_READ_FROM_MEMORY _particleVector includes only local particles, and in other case all (global) particles
	if(RRM_READ_FROM_MEMORY == _nReadMethod)
	{
		for(auto&& mol : _particleVector)
		{
			uint32_t cid_zb = mol.componentid();
			_density.at(cid_zb+1).numMolecules.local++;
		}
		// reduce
		uint16_t numComponents = global_simulation->getEnsemble()->getComponents()->size();
		domainDecomp->collCommInit(numComponents);
		for(uint16_t cid_ub=1; cid_ub<numComponents; ++cid_ub)
			domainDecomp->collCommAppendUnsLong(_density.at(cid_ub).numMolecules.local);

		domainDecomp->collCommAllreduceSum();

		for(uint16_t cid_ub=1; cid_ub<numComponents; ++cid_ub)
			_density.at(cid_ub).numMolecules.global = domainDecomp->collCommGetUnsLong();

		domainDecomp->collCommFinalize();
	}
	else
	{
		for(auto&& mol : _particleVector)
		{
			uint32_t cid_zb = mol.componentid();
			_density.at(cid_zb+1).numMolecules.global++;
		}
	}

	// sum up total number of particles, calc partial densities
	_density.at(0).numMolecules.global = 0;
	for(auto&& cid : _density)
	{
		_density.at(0).numMolecules.global += cid.numMolecules.global;
		cid.density = cid.numMolecules.global * dInvVolume;
	}
	_density.at(0).density = _density.at(0).numMolecules.global * dInvVolume;
}

void Reservoir::changeComponentID(Molecule& mol, const uint32_t& cid)
{
	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
	Component* compNew = &(ptrComps->at(_vecChangeCompIDs.at(cid) ) );
	mol.setComponent(compNew);
}

bool Reservoir::isRelevant(DomainDecompBase* domainDecomp, Domain* domain, Molecule& mol)
{
	double y = mol.r(1);
	uint32_t nBinIndex = floor(y / _dBinWidth);
	double dOffset;
	switch(_parent->getMovingDirection() )
	{
	case MD_LEFT_TO_RIGHT:
		dOffset = nBinIndex*_dBinWidth;
		break;
	case MD_RIGHT_TO_LEFT:
		dOffset = nBinIndex*_dBinWidth + (domain->getGlobalLength(1) - _dBinWidth);
		break;
	}
	return domainDecomp->procOwnsPos(mol.r(0), y-dOffset, mol.r(2), domain);
}

// queue methods
uint32_t Reservoir::getActualBinIndex() {return _binQueue->getActualBinIndex();}
uint64_t Reservoir::getNumMoleculesLocal() {return _binQueue->getNumParticles();}
uint32_t Reservoir::getNumBins() {return _binQueue->getNumBins();}
std::vector<Molecule>& Reservoir::getParticlesActualBin() {return _binQueue->getParticlesActualBin();}
bool Reservoir::nextBin(uint64_t& nMaxID) {bool bSuccess = _binQueue->next(); nMaxID += _density.at(0).numMolecules.global; return bSuccess;}
uint64_t Reservoir::getMaxMoleculeID() {return _binQueue->getMaxID();}
bool Reservoir::activateBin(uint32_t nBinIndex){return _binQueue->activateBin(nBinIndex);}
void Reservoir::clearBinQueue() {_binQueue->clear();}
void Reservoir::printBinQueueInfo()
{
	global_log->debug() << "_binQueue->getActualBinIndex()=" << _binQueue->getActualBinIndex() << endl;
	global_log->debug() << "_binQueue->getNumBins()=" << _binQueue->getNumBins() << endl;
	global_log->debug() << "_binQueue->getRoundCount()=" << _binQueue->getRoundCount() << endl;
	global_log->debug() << "_binQueue->getNumParticles()=" << _binQueue->getNumParticles() << endl;
	global_log->debug() << "_binQueue->getMaxID()=" << _binQueue->getMaxID() << endl;
}
