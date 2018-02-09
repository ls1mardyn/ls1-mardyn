/*
 * ParticleInsertion.cpp
 *
 *  Created on: 12.01.2018
 *      Author: mheinen
 */

#include "ParticleInsertion.h"
#include "utils/Random.h"
#include "utils/Region.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Component.h"
#include "molecules/Quaternion.h"
#include "NEMD/DensityControl.h"

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <vector>
#include <array>
#include <limits>
#include <algorithm>

template<typename T1, typename T2>
void select_rnd_elements(std::list<T1>& mylist, std::vector<T1>& myvec, T2 numSelect)
{
	myvec.clear();
	if(numSelect < 1)
		return;
	uint64_t numElements = mylist.size();
	uint64_t numElementsSub = numElements / numSelect;
	T2 numResidual = numElements % numSelect;

	std::vector<std::vector<T1> > mat;
	mat.resize(numSelect);
	for(T2 ci=0; ci<numSelect; ++ci)
	{
		if(ci<numResidual)
			mat.at(ci).resize(numElementsSub+1);
		else
			mat.at(ci).resize(numElementsSub);

		for(auto&& eli : mat.at(ci))
		{
			eli = mylist.front();
			mylist.pop_front();
		}

//		std::srand(std::time(nullptr));
		int rnd = rand()%numElementsSub;
		myvec.push_back(mat.at(ci).at(rnd) );
	}
}

// class ParticleManipDirector
ParticleManipDirector::ParticleManipDirector(dec::ControlRegion* region)
	:
	_region(region),
	_deleter(nullptr),
	_manipulator(nullptr)
{
	uint32_t numComps = global_simulation->getEnsemble()->getComponents()->size()+1;
	_inserter.resize(numComps);
//	_inserter.at(0) = nullptr;
	for(auto&& ins:_inserter)
		ins = new BubbleMethod(this, BMT_INSERTER);

	// deleter
	_deleter = new ParticleDeleter(this);

	// changer
	_changer.resize(numComps);
//	_changer.at(0) = nullptr;
	for(auto&& cit:_changer)
		cit = new BubbleMethod(this, BMT_CHANGER);

	// template methods
	DensityControl* ptrDC = dynamic_cast<DensityControl*>(region->GetParent() );
	_preForceAction = ptrDC->getPreForceAction();
	_postForceAction = ptrDC->getPostForceAction();
}

void ParticleManipDirector::readXML(XMLfileUnits& xmlconfig)
{
	string oldpath = xmlconfig.getcurrentnodepath();

	// set target variables for all components
	if(xmlconfig.changecurrentnode("targets")) {
		uint8_t numNodes = 0;
		XMLfile::Query query = xmlconfig.query("target");
		numNodes = query.card();
		global_log->info() << "Number of targets: " << (uint32_t)numNodes << endl;
		if(numNodes < 1) {
			global_log->error() << "No targets defined in XML-config file. Program exit ..." << endl;
			Simulation::exit(-1);
		}
		XMLfile::Query::const_iterator nodeIter;
		for( nodeIter = query.begin(); nodeIter; nodeIter++ ) {
			xmlconfig.changecurrentnode(nodeIter);
			int32_t cid = -1;
			xmlconfig.getNodeValue("@cid", cid);
			if(-1 == cid) {
				global_log->error() << "Missing attribute 'cid' in target. Program exit ..." << endl;
				Simulation::exit(-1);
			}
			else {
				// inserter
				if(xmlconfig.changecurrentnode("insertion") ) {
					_inserter.at(cid)->readXML(xmlconfig);
					_inserter.at(cid)->setTargetCompID(cid);
					xmlconfig.changecurrentnode("..");
				}

				// changer
				if(xmlconfig.changecurrentnode("change") ) {
					_changer.at(cid)->readXML(xmlconfig);
					_changer.at(cid)->setTargetCompID(cid);
					xmlconfig.changecurrentnode("..");
				}
			}
		}
	}
	xmlconfig.changecurrentnode(oldpath);
}

void ParticleManipDirector::localValuesReseted(Simulation* simulation)
{
	if(nullptr != _manipulator)
		_manipulator->Reset(simulation);
}

void ParticleManipDirector::globalValuesCalculated(Simulation* simulation)
{
	if(nullptr != _manipulator)
	{
		if(_manipulator->getManipState() == PMS_BUSY)
			return;
	}

	// global density > target value
	if(true == _region->globalTargetDensityExeeded(0) )
	{
		_manipulator = _deleter;
		global_log->info() << "DELETER activated" << endl;
	}
	else if(true == _region->globalTargetDensityUndershot(0) )
	{
		bool bValidIDs = this->setNextInsertIDs();
		if(true == bValidIDs)
		{
			_manipulator = _inserter.at(_nextChangeIDs.from);
//			_manipulator = nullptr;
			global_log->info() << "INSERTER activated" << endl;
		}
		else
			_manipulator = nullptr;
	}
	else
	{
		// check if particle identities have to be changed
		if(false == _region->globalCompositionBalanced() )
		{
			bool bValidIDs = this->setNextChangeIDs();
			if(true == bValidIDs)
			{
				_manipulator = _changer.at(_nextChangeIDs.from);
				global_log->info() << "CHANGER activated" << endl;
			}
			else
				_manipulator = nullptr;
		}
		else
			_manipulator = nullptr;
	}
	if(nullptr != _manipulator)
		_manipulator->PrepareParticleManipulation(simulation);
}

void ParticleManipDirector::ManipulateParticles(Simulation* simulation, Molecule* mol)
{
	if(nullptr != _manipulator)
		_manipulator->ManipulateParticles(simulation, mol);
}

void ParticleManipDirector::ManipulateParticleForces(Simulation* simulation, Molecule* mol)
{
	if(nullptr != _manipulator)
		_manipulator->ManipulateParticleForces(simulation, mol);
}

void ParticleManipDirector::FinalizeParticleManipulation(Simulation* simulation, MainLoopAction* action)
{
	if(nullptr != _manipulator)
	{
		if(action == _postForceAction)
			_manipulator->FinalizeParticleManipulation(simulation);
		else if(action == _preForceAction)
			_manipulator->FinalizeParticleManipulation_preForce(simulation);
	}
}

std::vector<dec::CompVarsStruct> ParticleManipDirector::getCompVars()
{
	return _region->getCompVars();
}

double ParticleManipDirector::GetLowerCorner(uint32_t nDim)
{
	return _region->GetLowerCorner(nDim);
}

double ParticleManipDirector::GetWidth(uint32_t nDim)
{
	return _region->GetWidth(nDim);
}

std::list<uint64_t> ParticleManipDirector::GetLocalParticleIDs(const uint32_t& nCompID)
{
	return _region->GetLocalParticleIDs(nCompID);
}

int64_t ParticleManipDirector::getLocalNumMoleculesSpread(uint32_t nCompID)
{
	return _region->getLocalNumMoleculesSpread(nCompID);
}

bool ParticleManipDirector::setNextChangeIDs()
{
	bool bRet = true;
	uint32_t numComps = global_simulation->getEnsemble()->getComponents()->size()+1;
	uint32_t cidMinSpread, cidMaxSpread;
	_region->getGlobalMinMaxNumMoleculesSpreadCompIDs(cidMinSpread, cidMaxSpread);
	_nextChangeIDs.from = cidMaxSpread;
	_nextChangeIDs.to   = cidMinSpread;
	bRet = bRet && (_nextChangeIDs.from > 0) && (_nextChangeIDs.from < numComps);
	bRet = bRet && (_nextChangeIDs.to   > 0) && (_nextChangeIDs.to   < numComps);
	return bRet;
}

bool ParticleManipDirector::setNextInsertIDs()
{
	bool bRet = true;
	uint32_t numComps = global_simulation->getEnsemble()->getComponents()->size()+1;
	uint32_t cidMinSpread, cidMaxSpread;
	_region->getGlobalMinMaxNumMoleculesSpreadCompIDs(cidMinSpread, cidMaxSpread);
	_nextChangeIDs.from = _nextChangeIDs.to = cidMinSpread;
	bRet = bRet && (_nextChangeIDs.from > 0) && (_nextChangeIDs.from < numComps);
	bRet = bRet && (_nextChangeIDs.to   > 0) && (_nextChangeIDs.to   < numComps);
	return bRet;
}

// class ParticleDeleter
ParticleDeleter::ParticleDeleter(ParticleManipDirector* director)
	:
	ParticleManipulator(director)
{
	uint32_t numComps = global_simulation->getEnsemble()->getComponents()->size()+1;
	_deletionLists.resize(numComps);
}

void ParticleDeleter::PrepareParticleManipulation(Simulation* simulation)
{
	_nManipState = PMS_BUSY;
	this->CreateDeletionLists(_director->getCompVars() );
}

void ParticleDeleter::ManipulateParticles(Simulation* simulation, Molecule* mol)
{
	ParticleContainer* particleCont = simulation->getMolecules();

	uint64_t mid = mol->id();
	uint64_t cid = mol->componentid()+1;
	bool bDeleteMolecule = false;
	for(auto did:_deletionLists.at(cid) )
	{
		if(did == mid)
			bDeleteMolecule = true;
	}
	if(true == bDeleteMolecule)
		particleCont->deleteMolecule(*mol, true);
}

void ParticleDeleter::CreateDeletionLists(std::vector<dec::CompVarsStruct> compVars)
{
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();
	uint8_t numComps = compVars.size();
	std::vector<CommVar<uint64_t> > numDel;
	numDel.resize(numComps);
	if(compVars.at(0).numMolecules.spread.global > 0)
		numDel.at(0).global = compVars.at(0).numMolecules.spread.global;
	else
	{
		for(uint8_t cid=0; cid<numComps; ++cid)
			_deletionLists.at(cid).clear();
		return;
	}
	std::vector<CommVar<int64_t> > positiveSpread;
	positiveSpread.resize(numComps);
	CommVar<int64_t> positiveSpreadSumOverComp;
	positiveSpreadSumOverComp.local  = 0;
	positiveSpreadSumOverComp.global = 0;

	for(uint8_t cid=1; cid<numComps; ++cid)
	{
		CommVar<int64_t> spread;
		// sum over components
		spread.global = compVars.at(cid).numMolecules.spread.global;
		if(spread.global > 0)
			positiveSpreadSumOverComp.global += spread.global;
		// sum over processes
		spread.local = compVars.at(cid).numMolecules.spread.local;
		if(spread.local > 0)
			positiveSpread.at(cid).local = spread.local;
		else
			positiveSpread.at(cid).local = 0;

	#ifdef ENABLE_MPI
		MPI_Allreduce( &positiveSpread.at(cid).local, &positiveSpread.at(cid).global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	#else
		positiveSpread.at(cid).global = positiveSpread.at(cid).local;
	#endif
	}

	std::vector<CommVar<double> > dInvPositiveSpread;
	dInvPositiveSpread.resize(numComps);
	CommVar<double> dInvPositiveSpreadSumOverComp;
	dInvPositiveSpreadSumOverComp.global = 1./( (double) (positiveSpreadSumOverComp.global) );

	for(uint8_t cid=1; cid<numComps; ++cid)
	{
		CommVar<int64_t> spread;
		// global
		spread.global = compVars.at(cid).numMolecules.spread.global;
		if(spread.global > 0)
			numDel.at(cid).global = round(spread.global * dInvPositiveSpreadSumOverComp.global * numDel.at(0).global);
		else
			numDel.at(cid).global = 0;
		// local
		dInvPositiveSpread.at(cid).global = 1./( (double) (positiveSpread.at(cid).global) );
		spread.local = compVars.at(cid).numMolecules.spread.local;
		cout << "Rank[" << domainDecomp.getRank() << "]: spread.local["<<(uint32_t)(cid)<<"] = " << spread.local << endl;
		if(spread.local > 0)
		{
			numDel.at(cid).local = round(spread.local * dInvPositiveSpread.at(cid).global * numDel.at(cid).global);
			cout << "Rank[" << domainDecomp.getRank() << "]: numDel.at("<<(uint32_t)(cid)<<").local = " << numDel.at(cid).local << endl;
			select_rnd_elements(compVars.at(cid).particleIDs, _deletionLists.at(cid), numDel.at(cid).local);
		}
		else
			numDel.at(cid).local = 0;

		//DEBUG
		cout << "Rank[" << domainDecomp.getRank() << "]compVars.at("<<(uint32_t)(cid)<<").deletionList:";
		for(auto mid:_deletionLists.at(cid) )
		{
			cout << " " << mid;
		}
		cout << endl;
		//DEBUG
	}
}

ParticleInsertion::ParticleInsertion(ParticleManipDirector* director, uint32_t state)
	:
	ParticleManipulator(director),
	_nState(state),
	_nTargetCompID(0)
{
}

BubbleMethod::BubbleMethod(ParticleManipDirector* director, uint32_t nType)
	:
	ParticleInsertion(director, BMS_IDLE),
	_rnd(nullptr),
	_selector(nullptr),
	_insRank(0),
	_maxID(100000000),
	_selectedMoleculeID(0),
	_nType(nType),
	_bVisual(false)
{
	// Random number generator
	_rnd = new Random(1000);

	// Particle selector
	_selector = new CompDependSelector(this);

	// init pbc vector
	Domain* domain = global_simulation->getDomain();
	double box[3];
	for(int16_t d=0; d<3; ++d)
		box[d] = domain->getGlobalLength(d);
	_pbc.initial.resize(27);
	int index = 0;
	for(int16_t xi=-1; xi<2; ++xi)
		for(int16_t yi=-1; yi<2; ++yi)
			for(int16_t zi=-1; zi<2; ++zi)
			{
				_pbc.initial.at(index).at(0) = xi*box[0];
				_pbc.initial.at(index).at(1) = yi*box[1];
				_pbc.initial.at(index).at(2) = zi*box[2];
				index++;
//				cout << index << ":("<<xi*box[0]<<","<<yi*box[1]<<","<<zi*box[2]<<")" << endl;
			}
}

BubbleMethod::~BubbleMethod()
{
	delete _rnd;
//	delete _selector; // TODO: not possible (protected)
}

void BubbleMethod::readXML(XMLfileUnits& xmlconfig)
{
	// visual
	{
		std::string strVisualUpper = "";
		bool bRet = xmlconfig.getNodeValue("@visual", strVisualUpper);
		transform(strVisualUpper.begin(), strVisualUpper.end(), strVisualUpper.begin(), ::toupper);
		_bVisual = bRet && ("YES" == strVisualUpper);
		if(true == _bVisual)
		{
			bRet = bRet && xmlconfig.getNodeValue("visIDs/selected", _visIDs.selected);
			bRet = bRet && xmlconfig.getNodeValue("visIDs/bubble", _visIDs.bubble);
			bRet = bRet && xmlconfig.getNodeValue("visIDs/force", _visIDs.force);
		}
		if(false == bRet)
		{
			global_log->error() << "BubbleMethod::readXML(): Something is wrong with visIDs! Program exit ..." << endl;
			Simulation::exit(-1);
		}
	}

	xmlconfig.getNodeValue("maxforce", _bubble.force.maxVal);
	xmlconfig.getNodeValue("forceradius", _bubble.forceRadius.target.global);
	_bubble.forceRadius.target.local = _bubble.forceRadius.actual.local = _bubble.forceRadius.target.global;
	_bubble.forceRadiusSquared.target.global = _bubble.forceRadiusSquared.target.global*_bubble.forceRadiusSquared.target.global;
	xmlconfig.getNodeValue("bubbleradius", _bubble.radius.target.global);
	_bubble.radius.target.local = _bubble.radius.actual.local = _bubble.radius.target.global;
	_bubble.radiusSquared.target.global = _bubble.radius.target.global*_bubble.radius.target.global;
	xmlconfig.getNodeValue("velocity/mean", _bubble.velocity.mean);
	xmlconfig.getNodeValue("velocity/max", _bubble.velocity.max);
	_bubble.velocity.maxSquared = _bubble.velocity.max * _bubble.velocity.max;

	if(this->checkBubbleRadius() == false)
	{
		global_log->error() << "BubbleMethod::readXML(): Bubble radius to large! Must be smaller than (0.5*boxLength - maxOuterMoleculeRadius)! Program exit ..." << endl;
		Simulation::exit(-1);
	}

	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();
	int ownRank = domainDecomp.getRank();
//	// DEBUG
//	cout << "rank[" << ownRank << "]: _bubble.radius.actual.global=" << _bubble.radius.actual.global << endl;
//	cout << "rank[" << ownRank << "]: _bubble.radius.target.global=" << _bubble.radius.target.global << endl;
//	cout << "rank[" << ownRank << "]: _bubble.radius.actual.local=" << _bubble.radius.actual.local << endl;
//	cout << "rank[" << ownRank << "]: _bubble.radius.target.local=" << _bubble.radius.target.local << endl;
//	// DEBUG

	// insertion molecules
	{
		string oldpath = xmlconfig.getcurrentnodepath();
		if(xmlconfig.changecurrentnode("molecules")) {
			uint8_t numNodes = 0;
			XMLfile::Query query = xmlconfig.query("molecule");
			numNodes = query.card();
			global_log->info() << "Number of insertion molecules: " << (uint32_t)numNodes << endl;
			if(numNodes < 1) {
				global_log->error() << "BubbleMethod::readXML(): No insertion molecules defined in XML-config file. Program exit ..." << endl;
				Simulation::exit(-1);
			}
			XMLfile::Query::const_iterator nodeIter;
			for( nodeIter = query.begin(); nodeIter; nodeIter++ ) {
				xmlconfig.changecurrentnode(nodeIter);
				double x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz;
				x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz = 0.0;
				xmlconfig.getNodeValue("x", x);
				xmlconfig.getNodeValue("y", y);
				xmlconfig.getNodeValue("z", z);
				xmlconfig.getNodeValue("vx", vx);
				xmlconfig.getNodeValue("vy", vy);
				xmlconfig.getNodeValue("vz", vz);
				xmlconfig.getNodeValue("q0", q0);
				xmlconfig.getNodeValue("q1", q1);
				xmlconfig.getNodeValue("q2", q2);
				xmlconfig.getNodeValue("q3", q3);
				xmlconfig.getNodeValue("Dx", Dx);
				xmlconfig.getNodeValue("Dy", Dy);
				xmlconfig.getNodeValue("Dz", Dz);
				Molecule mol(0, nullptr, x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz);
				_insertMolecules.initial.push_back(mol);
				_insertMolecules.actual.push_back(mol);
			}
		}
		xmlconfig.changecurrentnode(oldpath);
	}  // insertion molecules
}

bool BubbleMethod::checkBubbleRadius()
{
	std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
	double maxOuterMoleculeRadius = 0.0;
	for(auto cit:*ptrComps)
	{
		double tmp = cit.getOuterMoleculeRadiusLJ();
		if(tmp > maxOuterMoleculeRadius)
			maxOuterMoleculeRadius = tmp;
	}
	Domain* domain = global_simulation->getDomain();
	bool bValidRadius = true;
	for(uint8_t d=0; d<3; ++d)
		bValidRadius = (bValidRadius && (_bubble.radius.target.global < (0.5*domain->getGlobalLength(d)-maxOuterMoleculeRadius) ) );
	return bValidRadius;
}

void BubbleMethod::Reset(Simulation* simulation)
{
	// reset local values
	_numManipulatedParticles.local = 0;
	_bubble.radius.actual.local = _bubble.radius.target.global;

	ParticleContainer* particleCont = simulation->getMolecules();
	ParticleIterator pit;
	for( pit  = particleCont->iteratorBegin();
			pit != particleCont->iteratorEnd();
		 ++pit )
	{
		this->resetBubbleMoleculeComponent(&(*pit) );
	}
	// reset bubble molecules vector
	_bubbleMolecules.clear();
}

void BubbleMethod::PrepareParticleManipulation(Simulation* simulation)
{
	_nManipState = PMS_BUSY;

	// reset local values
	this->Reset(simulation);

	ChangeVar<uint32_t> nextChangeIDs = _director->getNextChangeIDs();
	_nextChangeIDs.from = nextChangeIDs.from;
	_nextChangeIDs.to = nextChangeIDs.to;
	this->selectParticle(simulation);
	if(BMT_INSERTER == _nType)
		this->initInsertionMolecules(simulation);
}

void BubbleMethod::ManipulateParticles(Simulation* simulation, Molecule* mol)
{
	this->FreezeSelectedMolecule(mol);
}

void BubbleMethod::ManipulateParticleForces(Simulation* simulation, Molecule* mol)
{
	this->GrowBubble(simulation, mol);
}

void BubbleMethod::FreezeSelectedMolecule(Molecule* mol)
{
	if(_nState != BMS_GROWING_BUBBLE)
		return;

	if(mol->id() == _selectedMoleculeID)
	{
		for(uint8_t d=0; d<3; ++d)
			mol->setr(d, _selectedMoleculeInitPos.at(d) );
	}
}

void BubbleMethod::resetBubbleMoleculeComponent(Molecule* mol)
{
	if(_nState != BMS_GROWING_BUBBLE)
		return;

	for(auto mit:_bubbleMolecules)
	{
		if(mol->id() == mit.mid)
		{
			Component* comp = global_simulation->getEnsemble()->getComponent(mit.cid);
			mol->setComponent(comp);
		}
	}
}

double BubbleMethod::calcMinSquaredDistance(Molecule* mol, const double& radius)
{
	double dOuterMoleculeRadiusLJ = mol->component()->getOuterMoleculeRadiusLJ();
	double dist2_min = (radius+dOuterMoleculeRadiusLJ)*(radius+dOuterMoleculeRadiusLJ);
	return dist2_min;
}

bool BubbleMethod::outerMoleculeRadiusCutsBubbleRadius(Simulation* simulation, Molecule* mol, const double& dist2_min, double& dist2, double* distVec)
{
	DomainDecompBase domainDecomp = simulation->domainDecomposition();
	int ownRank = domainDecomp.getRank();
	Domain* domain = simulation->getDomain();
	double box[3];
	for(uint8_t d=0; d<3; ++d)
		box[d] = domain->getGlobalLength(d);

	bool bCutsBubbleRadius = false;
	for(auto pbc:_pbc.initial)
	{
		dist2 = 0.0;
		for(uint8_t d=0; d<3; ++d)
		{
			distVec[d] = mol->r(d) - (_selectedMoleculeInitPos.at(d)+pbc.at(d) );
			dist2 += distVec[d]*distVec[d];
		}
		if(dist2 < dist2_min)
		{
			bCutsBubbleRadius = true;
			break;
		}
	}
	return bCutsBubbleRadius;
}

void BubbleMethod::updateActualBubbleRadiusLocal(Molecule* mol, const double& dist)
{
	double tmp = dist - mol->component()->getOuterMoleculeRadiusLJ();
	if(tmp < _bubble.radius.actual.local)
		_bubble.radius.actual.local = tmp;
}

void BubbleMethod::updateForceOnBubbleCuttingMolecule(Molecule* mol, const double& dist2_min, const double& dist2, const double& dist, double* distVec)
{
	double v2 = mol->v2();
	if(v2 >= _bubble.velocity.maxSquared)
		return;
	double invDist = 1./dist;
	double Fadd[3];
	for(uint8_t d=0; d<3; ++d)
		Fadd[d] = distVec[d] * _bubble.force.maxVal * (1-dist2/dist2_min);
	mol->Fadd(Fadd);
}

void BubbleMethod::GrowBubble(Simulation* simulation, Molecule* mol)
{
	DomainDecompBase domainDecomp = simulation->domainDecomposition();
	int ownRank = domainDecomp.getRank();
//	cout << "rank[" << ownRank << "]: GrowBubble(), state:" << _nState << endl;

	if(_nState != BMS_GROWING_BUBBLE)
	{
//		cout << "postEventNewTimestepAction(), state:" << _nState << endl;
		return;
	}

	/* DEBUG -->
	double Fabs = 0.0;
	for(uint8_t dim=0; dim<3; ++dim)
		Fabs += mol->F(dim)*mol->F(dim);
	Fabs = sqrt(Fabs);
	if(Fabs > 350)
	{
		cout << "rank[" << ownRank << "]: id=" << mol->id() << ", Fabs=" << Fabs << endl;

		Component* comp6 = global_simulation->getEnsemble()->getComponent(5);  // max_force_exceed
		ID_pair ids;
		ids.mid = mol->id();
		ids.cid = mol->componentid();
		_bubbleMolecules.push_back(ids);
		mol->setComponent(comp6);
	}
	// <-- DEBUG */

	// skip selected molecule, but update its position
	if(mol->id() == _selectedMoleculeID)
	{
		cout << "rank[" << ownRank << "] BubbleMethod::GrowBubble: _selectedMoleculeID="<<_selectedMoleculeID<< endl;
		cout << *mol << endl;
		// visual
		if(true == _bVisual)
		{
			Component* comp = global_simulation->getEnsemble()->getComponent(_visIDs.selected-1);  // selected
			ID_pair ids;
			ids.mid = mol->id();
			ids.cid = mol->componentid();
			_bubbleMolecules.push_back(ids);
			mol->setComponent(comp);
		}

		return;
	}

	// FORCE RADIUS -->
	bool bCutsForceRadius;
	{
		double distVec[3] = {0.0, 0.0, 0.0};
		double dist2 = 0.0;
		double dist2_min = calcMinSquaredDistance(mol, _bubble.forceRadius.target.global);
		bCutsForceRadius = this->outerMoleculeRadiusCutsBubbleRadius(simulation, mol, dist2_min, dist2, distVec);

		if(false == bCutsForceRadius)
			return;

		double dist = sqrt(dist2);
		this->updateForceOnBubbleCuttingMolecule(mol, dist2_min, dist2, dist, distVec);
		_numManipulatedParticles.local++;
	}  // <-- FORCE RADIUS

	// BUBBLE RADIUS -->
	{
		double distVec[3] = {0.0, 0.0, 0.0};
		double dist2 = 0.0;
		double dist2_min = calcMinSquaredDistance(mol, _bubble.radius.target.global);
		bool bCutsBubbleRadius = this->outerMoleculeRadiusCutsBubbleRadius(simulation, mol, dist2_min, dist2, distVec);

		/* DEBUG -->
		if(mol->id() == 100000382)
		{
			cout << "rank[" << ownRank << "]: GrowBubble(), ID=" << mol->id() << ", (x,y,z): " << mol->r(0) << "," << mol->r(1) << "," << mol->r(2) << endl;
			cout << "rank[" << ownRank << "]: GrowBubble(), dist2=" << dist2 << ", dist=" << sqrt(dist2) << endl;
		}
		// <-- DEBUG */

		if(true == bCutsForceRadius)
		{
			uint32_t cid = _visIDs.force-1;
			if(true == bCutsBubbleRadius)
			{
				cid = _visIDs.bubble-1;
				double dist = sqrt(dist2);
				this->updateActualBubbleRadiusLocal(mol, dist);
			}
			Component* comp = global_simulation->getEnsemble()->getComponent(cid);
			ID_pair ids;
			ids.mid = mol->id();
			ids.cid = mol->componentid();
			_bubbleMolecules.push_back(ids);
			mol->setComponent(comp);
		}
	}  // <-- BUBBLE RADIUS
}

void BubbleMethod::ChangeIdentity(Simulation* simulation, Molecule* mol)
{
	// --> CHANGE_IDENTITY
	{
		std::vector<Component>* ptrComps = simulation->getEnsemble()->getComponents();
		Component* compOld = mol->component();
		uint32_t newID = _nextChangeIDs.to-1;
		if(newID > ptrComps->size() )
		{
			global_log->error() << "Component ID (1st Component: ID=1): " << newID+1 << " does NOT exist! Program exit ..." << endl;
			Simulation::exit(-1);
		}
		Component* compNew = &(ptrComps->at(newID) );

		uint8_t numRotDOF_old = compOld->getRotationalDegreesOfFreedom();
		uint8_t numRotDOF_new = compNew->getRotationalDegreesOfFreedom();
		double dUkinOld = mol->U_kin();
		double dUkinPerDOF = dUkinOld / (3 + numRotDOF_old);

		// rotation
		double U_rot = mol->U_rot();
		global_log->info() << "U_rot_old = " << U_rot << endl;

		double L[3];
		double Ipa[3];

		Ipa[0] = compNew->I11();
		Ipa[1] = compNew->I22();
		Ipa[2] = compNew->I33();

		for(uint8_t dim=0; dim<3; ++dim)
		{
			L[dim] = sqrt(dUkinPerDOF * 2. * Ipa[dim] );
			mol->setD(dim, L[dim] );
		}

#ifndef NDEBUG
		cout << "L[0] = " << L[0] << endl;
		cout << "L[1] = " << L[1] << endl;
		cout << "L[2] = " << L[2] << endl;
#endif
		Quaternion q(1., 0., 0., 0.);
		mol->setq(q);

		mol->setComponent(compNew);
//		mol->clearFM();  // <-- necessary?
#ifndef NDEBUG
		cout << "Changed cid of molecule " << mol->id() << " from: " << _nextChangeIDs.from << " to: " << mol->componentid()+1 << endl;
#endif

		// update transl. kin. energy
		double dScaleFactorTrans = sqrt(6*dUkinPerDOF/compNew->m()/mol->v2() );
		mol->scale_v(dScaleFactorTrans);

		U_rot = mol->U_rot();
		global_log->info() << "U_rot_new = " << U_rot << endl;

		/*
		//connection to MettDeamon
		if(NULL != _mettDeamon)
			_mettDeamon->IncrementChangedMoleculesLocal();
		*/
	}
	// <-- CHANGE_IDENTITY
}

void BubbleMethod::FinalizeParticleManipulation(Simulation* simulation)
{
	DomainDecompBase domainDecomp = simulation->domainDecomposition();
	int ownRank = domainDecomp.getRank();
	ParticleContainer* particleCont = global_simulation->getMolecules();

#ifdef ENABLE_MPI
	MPI_Allreduce( &_numManipulatedParticles.local, &_numManipulatedParticles.global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce( &_bubble.radius.actual.local, &_bubble.radius.actual.global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
	_numManipulatedParticles.global = _numManipulatedParticles.local;
	_bubble.radius.actual.global = _bubble.radius.actual.local;
#endif

//	cout << "rank[" << ownRank << "]: _numManipulatedParticles.global=" << _numManipulatedParticles.global << endl;
	cout << "rank[" << ownRank << "]: _bubble.radius.actual.global=" << _bubble.radius.actual.global << endl;
//	cout << "rank[" << ownRank << "]: _bubble.radius.target.global=" << _bubble.radius.target.global << endl;
//	cout << "rank[" << ownRank << "]: BubbleMethod::FinalizeParticleManipulation(): _selectedMoleculeID="<<_selectedMoleculeID<< endl;
//	cout << "rank[" << ownRank << "]: BubbleMethod::FinalizeParticleManipulation(): _nState="<<_nState<< endl;

//	if( !(BMS_GROWING_BUBBLE == _nState && _numManipulatedParticles.global == 0) )

	bool bBubbleSizeReached = (_bubble.radius.actual.global >= _bubble.radius.target.global) && BMS_GROWING_BUBBLE == _nState;
//	cout << "rank[" << ownRank << "]: BubbleMethod::FinalizeParticleManipulation(): bBubbleSizeReached="<<bBubbleSizeReached<< ", targetID=" << this->getTargetCompID() << endl;
	if(false == bBubbleSizeReached)
		return;

	// change back components, TODO: maybe better location to do so, to avoid additional loop
	{
		ParticleContainer* particles = global_simulation->getMolecules();
		ParticleIterator pit;
		for( pit  = particles->iteratorBegin();
				pit != particles->iteratorEnd();
			 ++pit )
		{
			this->resetBubbleMoleculeComponent( &(*pit) );
		}
	}

	if(ownRank == _insRank)
	{
		if(BMT_CHANGER == this->getType() )
		{
			ParticleContainer* particles = global_simulation->getMolecules();
			ParticleIterator pit;
			for( pit  = particles->iteratorBegin();
					pit != particles->iteratorEnd();
				 ++pit )
			{
				if(pit->id() == _selectedMoleculeID) //componentid()+1 == 3)
				{
					this->ChangeIdentity(simulation, &(*pit) );
				}
			}
		}
		else if(BMT_INSERTER == this->getType() )
		{
			ParticleContainer* particles = global_simulation->getMolecules();
			ParticleIterator pit;
			for( pit  = particles->iteratorBegin();
					pit != particles->iteratorEnd();
				 ++pit )
			{
				this->resetBubbleMoleculeComponent( &(*pit) );
				if(pit->id() == _selectedMoleculeID) //componentid()+1 == 3)
				{
					particleCont->deleteMolecule(*pit, false);
					break;
				}
			}
	//		particleCont->addParticles(_insertMolecules, false);
			for(auto&& mol:_insertMolecules.actual)
			{
				cout << "rank[" << ownRank << "]: Adding particle..." << endl;
				cout << mol;
				particleCont->addParticle(mol, true, true, true);
				cout << "rank[" << ownRank << "]: ... added!" << endl;
			}
		}
	}
	// update maxID
#ifdef ENABLE_MPI
	MPI_Bcast( &_maxID, 1, MPI_UNSIGNED_LONG, _insRank, MPI_COMM_WORLD);
#endif
	_nState = BMS_IDLE;
	_nManipState = PMS_IDLE;
	return;

	/*
	 * change component for visualization
	 *
	 */
	Component* comp2 = global_simulation->getEnsemble()->getComponent(1);  // H2
	Component* comp6 = global_simulation->getEnsemble()->getComponent(5);  // H2 insert

	ParticleContainer* particles = global_simulation->getMolecules();
	ParticleIterator pit;
	double dBubbleRadius = _bubble.radius.target.global;
	double vm = _bubble.velocity.mean;
	for( pit  = particles->iteratorBegin();
			pit != particles->iteratorEnd();
		 ++pit )
	{
		if(pit->id() == _selectedMoleculeID) //componentid()+1 == 3)
		{
//			pit->setComponent(comp6);
			pit->setr(0, pit->r(0)-dBubbleRadius*0.5);
			pit->setv(0, -vm);
		}
//		else if(pit->componentid()+1 == 4)
//			pit->setComponent(comp2);
	}

	if(ownRank == _insRank)
	{
		double rx, ry, rz;
		rx = _selectedMoleculeInitPos.at(0);
		ry = _selectedMoleculeInitPos.at(1);
		rz = _selectedMoleculeInitPos.at(2);
		Molecule mol1(++_maxID, comp6, rx+dBubbleRadius*0.5, ry, rz, vm, 0, 0);
		Molecule mol2(++_maxID, comp6, rx, ry+dBubbleRadius*0.5, rz, 0,  vm, 0);
		Molecule mol3(++_maxID, comp6, rx, ry-dBubbleRadius*0.5, rz, 0, -vm, 0);
		Molecule mol4(++_maxID, comp6, rx, ry, rz+dBubbleRadius*0.5, 0, 0,  vm);
		Molecule mol5(++_maxID, comp6, rx, ry, rz-dBubbleRadius*0.5, 0, 0, -vm);

		cout << "rank[" << ownRank << "]: postLoopAction() adding particle at: " << rx << ", " << ry << ", " << rz << endl;
		particleCont->addParticle(mol1, true, true, true);
		particleCont->addParticle(mol2, true, true, true);
		particleCont->addParticle(mol3, true, true, true);
		particleCont->addParticle(mol4, true, true, true);
		particleCont->addParticle(mol5, true, true, true);
	}

	// update maxID
#ifdef ENABLE_MPI
	MPI_Bcast( &_maxID, 1, MPI_UNSIGNED_LONG, _insRank, MPI_COMM_WORLD);
#endif

	_nState = BMS_IDLE;
}

void BubbleMethod::FinalizeParticleManipulation_preForce(Simulation* simulation)
{
}

void BubbleMethod::selectParticle(Simulation* simulation)
{
	if(_nState != BMS_IDLE)
		return;

	DomainDecompBase* domainDecomp = &(simulation->domainDecomposition() );
	Molecule selectedMolecule;
	_insRank = _selector->selectParticle(simulation, selectedMolecule);
//	Component* comp4 = global_simulation->getEnsemble()->getComponent(3);  // selected
//	selectedMolecule.setComponent(comp4);
	uint64_t selectedMoleculeID = 0;
	double dPos[3] = {0.0, 0.0, 0.0};
//	cout << "BubbleMethod::selectParticle: selectedMolecule.id()=" << selectedMolecule.id() << endl;
	if(domainDecomp->getRank() == _insRank)
	{
//		if(nullptr == selectedMolecule)
//		{
//			global_log->error() << "BubbleMethod::selectParticle: Condition (nullptr == selectedMolecule) failed! Programm exit ..." << endl;
//			Simulation::exit(-1);
//		}
		selectedMoleculeID = selectedMolecule.id();
		for(uint8_t d=0; d<3; ++d)
			dPos[d] = selectedMolecule.r(d);
	}
#ifdef ENABLE_MPI
	MPI_Bcast( &selectedMoleculeID, 1, MPI_UNSIGNED_LONG, _insRank, MPI_COMM_WORLD);
	MPI_Bcast( dPos, 3, MPI_DOUBLE, _insRank, MPI_COMM_WORLD);
#endif
	// store selected molecule values
	_selectedMoleculeID = selectedMoleculeID;
	for(uint8_t d=0; d<3; ++d)
		_selectedMoleculeInitPos.at(d) = dPos[d];

	// inform director about selected particle
//	_director->
	_nState = BMS_GROWING_BUBBLE;
}

Quaternion BubbleMethod::createRandomQuaternion()
{
	double alpha_rad_half = M_PI*_rnd->rnd();
	std::array<double,3> n;
	for(uint16_t d=0; d<3; ++d)
		n.at(d) = _rnd->rnd();
	Quaternion q(alpha_rad_half, n);
//	cout << "q=("<< q.qw() << "," << q.qx() << "," << q.qy() << "," << q.qz() << ")" << endl;
	return q;
}

void BubbleMethod::initInsertionMolecules(Simulation* simulation)
{
	for(size_t mi=0; mi<_insertMolecules.actual.size(); ++mi)
	{
		for(uint8_t d=0; d<3; ++d)
			_insertMolecules.actual.at(mi).setr(d, _insertMolecules.initial.at(mi).r(d) );
	}

	std::vector<Component>* ptrComps = simulation->getEnsemble()->getComponents();
	Component* compIns = &(ptrComps->at(_nTargetCompID-1) );

	Quaternion q = this->createRandomQuaternion();
	Molecule selectedMolecule;
	for(uint16_t d=0; d<3; ++d)
		selectedMolecule.setr(d, _selectedMoleculeInitPos.at(d) );
	for(auto&& mol:_insertMolecules.actual)
	{
		// rotate insertion molecules around insertion position -->
		std::array<double, 3> pos;
		for(uint16_t d=0; d<3; ++d)
			pos.at(d) = mol.r(d);
		q.rotateInPlace(pos);
		for(uint16_t d=0; d<3; ++d)
			pos.at(d) += _selectedMoleculeInitPos.at(d);
		// <-- rotate insertion molecules around insertion position

		// set position and component
		mol.setid(++_maxID);
		mol.setComponent(compIns);
		for(uint16_t d=0; d<3; ++d)
			mol.setr(d, pos.at(d) );

		// set velocity, so that inserted pair of molecules moves apart: <-- O ... O -->
		double vi[3];
		double dist2 = selectedMolecule.dist2(mol, vi);
		double inv_dist = 1./sqrt(dist2);
		for(uint16_t d=0; d<3; ++d)
			mol.setv(d, vi[d]*inv_dist*_bubble.velocity.mean);
	}
	// DEBUG -->
	double dr_initial[3];
	double dr_actual[3];
	double dist2_initial = _insertMolecules.initial.at(0).dist2(_insertMolecules.initial.at(1), dr_initial);
	double dist2_actual = _insertMolecules.actual.at(0).dist2(_insertMolecules.actual.at(1), dr_actual);
	double dist2_diff = dist2_actual - dist2_initial;
	mardyn_assert(fabs(dist2_diff) <= 1e-15);
	// <-- DEBUG
}

// Selector info
double BubbleMethod::GetLowerCorner(uint32_t nDim)
{
	return _director->GetLowerCorner(nDim);
}

double BubbleMethod::GetWidth(uint32_t nDim)
{
	return _director->GetLowerCorner(nDim);
}

std::list<uint64_t> BubbleMethod::GetLocalParticleIDs(const uint32_t& nCompID)
{
	return _director->GetLocalParticleIDs(nCompID);
}

int64_t BubbleMethod::getLocalNumMoleculesSpread(uint32_t nCompID)
{
	return _director->getLocalNumMoleculesSpread(nCompID);
}

// class RandomSelector : public ParticleSelector
RandomSelector::RandomSelector(BubbleMethod* parent)
	:
	ParticleSelector(parent),
	_random(nullptr)
{
	_random = new Random;
}

void RandomSelector::generateRandomInsPos()
{
	for(uint8_t d=0; d<3; ++d) {
		float rnd = _random->rnd();
		_daInsertionPosition.at(d) = _insRegLowerCorner.at(d) + rnd * _insRegWidth.at(d);
	}
}

int RandomSelector::selectParticle(Simulation* simulation, Molecule& selectedMolecule)
{
	this->collectInfo();

	DomainDecompBase domainDecomp = simulation->domainDecomposition();
	int ownRank = domainDecomp.getRank();
	cout << "rank[" << ownRank << "]: selectParticle()" << endl;

	Domain* domain = global_simulation->getDomain();
	double bbMin[3];
	double bbMax[3];

	domainDecomp.getBoundingBoxMinMax(domain, bbMin, bbMax);

	this->generateRandomInsPos();

	ParticleContainer* particles = simulation->getMolecules();
	uint64_t numMols = particles->getNumberOfParticles();

	// find molecule with min distance to insertion position
	cout << "rank[" << ownRank << "]: find molecule with min distance to insertion position" << endl;
	double dist2_min = 1000;

	ParticleIterator pit;
	for( pit  = particles->iteratorBegin();
			pit != particles->iteratorEnd();
		 ++pit )
	{
		double dist2 = 0.0;
		double distVec[3];
		for(uint8_t d=0; d<3; ++d)
		{
			distVec[d] = pit->r(d) - _daInsertionPosition.at(d);
			dist2 += distVec[d]*distVec[d];
		}
		if(dist2 < dist2_min)
		{
			dist2_min = dist2;
			selectedMolecule = *pit;
		}
	}

	struct {
		double val;
		int   rank;
	} in, out;

	in.val = dist2_min;
	in.rank = ownRank;

#ifdef ENABLE_MPI
	MPI_Allreduce( &in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
#else
	out.rank = ownRank;
#endif
	cout << "rank["<<ownRank<<"]: dist2_min="<<dist2_min<<endl;

	return out.rank;
}

void RandomSelector::collectInfo()
{
	for(uint8_t d=0; d<3; ++d)
	{
		_insRegLowerCorner.at(d) = _parent->GetLowerCorner(d);
		_insRegWidth.at(d) = _parent->GetWidth(d);
	}
}


// class CompDependSelector : public ParticleSelector
CompDependSelector::CompDependSelector(BubbleMethod* parent)
	:
	ParticleSelector(parent)
{
}

int CompDependSelector::selectParticle(Simulation* simulation, Molecule& selectedMolecule)
{
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();
	int ownRank = domainDecomp.getRank();
	cout << "ownRank=" << ownRank << endl;
	ChangeVar<uint32_t> nextChangeIDs = _parent->getNextChangeIDs();
	std::list<uint64_t> particleIDsFrom = _parent->GetLocalParticleIDs(nextChangeIDs.from);
	std::list<uint64_t> particleIDsTo   = _parent->GetLocalParticleIDs(nextChangeIDs.to);
	int64_t spreadDiff = std::numeric_limits<int64_t>::min();
	if(particleIDsFrom.size() > 0)
	{
		int64_t spreadFrom = _parent->getLocalNumMoleculesSpread(nextChangeIDs.from);
		int64_t spreadTo   = _parent->getLocalNumMoleculesSpread(nextChangeIDs.to);
		if(nextChangeIDs.from != nextChangeIDs.to)  // change
			spreadDiff = spreadFrom - spreadTo;
		else
			spreadDiff = spreadFrom*-1;  // insert

		/* DEBUG -->
		cout << "Rank[" << ownRank << "]:  spreadFrom=" << spreadFrom
			<< ", spreadTo=" << spreadTo
			<< ", spreadDiff=" << spreadDiff << endl;
		// <-- DEBUG */
	}
	cout << "Rank[" << ownRank << "]:  spreadDiff=" << spreadDiff << endl;

	struct {
		int64_t val;
		int    rank;
	} in, out;

	in.val = spreadDiff;
	in.rank = ownRank;

#ifdef ENABLE_MPI
	MPI_Allreduce( &in, &out, 1, MPI_LONG_INT, MPI_MAXLOC, MPI_COMM_WORLD);
#else
	out.rank = ownRank;
#endif

	uint64_t selectedID = 0;
	if(ownRank == out.rank)
	{
		// DEBUG -->
		cout << "Rank[" << ownRank << "]: particleIDsFrom:" << endl;
		for(auto&& id:particleIDsFrom)
			cout << id << ", ";
		cout << endl;
		// <-- DEBUG
		std::vector<uint64_t> selectedIDs;
		if(particleIDsFrom.size() < 1)
		{
			global_log->error() << "CompDependSelector::selectParticle: Condition (particleIDsFrom.size() < 1) failed! Program exit ..." << endl;
			Simulation::exit(-1);
		}
		select_rnd_elements(particleIDsFrom, selectedIDs, 1);
		if(selectedIDs.size() != 1)
		{
			global_log->error() << "CompDependSelector::selectParticle: Condition (selectedIDs.size() != 1) failed! Program exit ..." << endl;
			Simulation::exit(-1);
		}
		selectedID = selectedIDs.at(0);
	}

	ParticleContainer* particles = simulation->getMolecules();
	ParticleIterator pit;
	for( pit  = particles->iteratorBegin();
			pit != particles->iteratorEnd();
		 ++pit )
	{
		uint64_t id = pit->id();
		if(selectedID == id)
		{
			selectedMolecule = *pit;
//			cout << "selectedMolecule=" << selectedMolecule << endl;
			break;
		}
	}
//	cout << "out.rank=" << out.rank << endl;
//	cout << "selectedMolecule.id()=" << selectedMolecule.id() << endl;
	return out.rank;
}

void CompDependSelector::collectInfo()
{

}
