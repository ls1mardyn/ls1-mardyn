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

#include <cstdlib>
#include <iostream>
#include <ctime>

ParticleInsertion::ParticleInsertion(CuboidRegion* parent, uint32_t state)
	:
	_parent(parent),
	_nState(state)
{
}

BubbleMethod::BubbleMethod(CuboidRegion* parent)
	:
	ParticleInsertion(parent, BMS_IDLE),
	_random(nullptr),
	_dBubbleForce(0.0),
	_dBubbleRadius(0.0),
	_dBubbleRadiusSquared(0.0),
	_dMoleculeOuterRadius(0.0),
	_dReplacement(0.0),
	_dMaxForce(0.0),
	_vm(1.0),
	_insRank(0),
	_numReplacements(0),
	_maxID(100000000),
	_selectedMoleculeID(0)
{
	_random = new Random();
}

BubbleMethod::~BubbleMethod()
{
	delete _random;
}

void BubbleMethod::readXML(XMLfileUnits& xmlconfig)
{
	xmlconfig.getNodeValue("replacement", _dReplacement);
	xmlconfig.getNodeValue("maxforce", _dMaxForce);
	xmlconfig.getNodeValue("molouterradius", _dMoleculeOuterRadius);
	xmlconfig.getNodeValue("bubbleradius", _dBubbleRadius);
	_dBubbleRadiusSquared = _dBubbleRadius*_dBubbleRadius;
	xmlconfig.getNodeValue("vm", _vm);
}

void BubbleMethod::preLoopAction()
{
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();
	int rank = domainDecomp.getRank();
	cout << "rank[" << rank << "]: preLoopAction(), state:" << _nState << endl;
	// Reset replacement counter
	_numReplacements = 0;

	if(_nState != BMS_SELECT_RANDOM_MOLECULE)
		return;

	Domain* domain = global_simulation->getDomain();
	double bbMin[3];
	double bbMax[3];
	double insPos[3];
	uint64_t insID;
	domainDecomp.getBoundingBoxMinMax(domain, bbMin, bbMax);

	for(uint8_t d=0; d<3; ++d)
	{
		float rnd = _random->rnd();
//		cout << "rnd=" << rnd << endl;
//		_daInsertionPosition.at(d) = rnd * (bbMax[d] - bbMin[d]);
		insPos[d] = _parent->GetLowerCorner(d) + rnd * _parent->GetWidth(d);
	}
	// DEBUG
	_daInsertionPosition.at(2) = global_simulation->getDomain()->getGlobalLength(2) * 0.5;
	// DEBUG

	ParticleContainer* particles = global_simulation->getMolecules();
	uint64_t numMols = particles->getNumberOfParticles();

	/*
	std::srand(std::time(nullptr)); // use current time as seed for random generator
	int rnd = std::rand() % numMols;
	ParticleIterator tM  = particles->iteratorBegin() += rnd;
	*/

	// find molecule with min distance to insertion position
	cout << "rank[" << rank << "]: find molecule with min distance to insertion position" << endl;
	double dist2_min = 1000;
	Molecule* selectedMolecule = nullptr;
	ParticleIterator pit;
	for( pit  = particles->iteratorBegin();
			pit != particles->iteratorEnd();
		 ++pit )
	{
		double dist2 = 0.0;
		double distVec[3];
		for(uint8_t d=0; d<3; ++d)
		{
			distVec[d] = pit->r(d) - insPos[d];
			dist2 += distVec[d]*distVec[d];
		}
		if(dist2 < dist2_min)
		{
			dist2_min = dist2;
			selectedMolecule = &(*pit);
		}
	}

	double aout;
	int ind;
	struct {
		double val;
		int   rank;
	} in, out;

	in.val = dist2_min;
	in.rank = rank;

	MPI_Allreduce( &in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
	_insRank = out.rank;

	cout << "rank["<<rank<<"]: dist2_min="<<dist2_min<<endl;
	if(rank == _insRank)
	{
		cout << "------------------------------------" << endl;
		cout << "dist2_min="<<out.val<<endl;
		cout << "rank_min="<<out.rank<<endl;
		cout << "------------------------------------" << endl;

		if(nullptr == selectedMolecule)
			return;

		Component* comp4 = global_simulation->getEnsemble()->getComponent(3);  // H2 mark selected molecule
		selectedMolecule->setComponent(comp4);

		insID = selectedMolecule->id();
		for(uint8_t d=0; d<3; ++d)
			insPos[d] = selectedMolecule->r(d);
		cout << "Selected molecule: id=" << _selectedMoleculeID << ", pos: " << selectedMolecule->r(0) << ", " << selectedMolecule->r(1) << ", " << selectedMolecule->r(2) << endl;
	}
	// Broadcast
	MPI_Bcast( insPos, 3, MPI_DOUBLE, out.rank, MPI_COMM_WORLD);
	MPI_Bcast( &insID, 1, MPI_UNSIGNED_LONG, out.rank, MPI_COMM_WORLD);

	// store data
	_selectedMoleculeID = insID;
	for(uint8_t d=0; d<3; ++d)
		_selectedMoleculePos.at(d) = insPos[d];

	cout << "rank["<<rank<<"]: _selectedMoleculeID="<<_selectedMoleculeID<<endl;

	_nState = BMS_GROWING_BUBBLE;
}

void BubbleMethod::insideLoopAction(Molecule* mol)
{
	if(_nState != BMS_GROWING_BUBBLE)
	{
//		cout << "insideLoopAction(), state:" << _nState << endl;
		return;
	}

	// skip selected molecule, but update its position
//	cout << "_selectedMoleculeID="<<_selectedMoleculeID<<",mol->id()="<<mol->id()<<endl;
	if(mol->id() == _selectedMoleculeID)
	{
//		for(uint8_t d=0; d<3; ++d)
//			_selectedMoleculePos.at(d) = mol->r(d);
//		cout << "insideLoopAction(): molecule: id=" << mol->id() << ", pos: " << mol->r(0) << ", " << mol->r(1) << ", " << mol->r(2) << endl;
		return;
	}
//	else
//		global_log->error() << "Selected molecule NOT found!" << endl;

	double dist2 = 0.0;
	double distVec[3];
	for(uint8_t d=0; d<3; ++d)
	{
		distVec[d] = mol->r(d) - _selectedMoleculePos.at(d);
		dist2 += distVec[d]*distVec[d];
	}

	double dist2_min = (_dBubbleRadius+_dMoleculeOuterRadius)*(_dBubbleRadius+_dMoleculeOuterRadius);
	Component* comp2 = global_simulation->getEnsemble()->getComponent(1);  // H2
	Component* comp5 = global_simulation->getEnsemble()->getComponent(4);  // H2 replace
	if(mol->componentid()+1 == 5)
		mol->setComponent(comp2);

	if(dist2 > dist2_min)
		return;
	else
		mol->setComponent(comp5);

	double dist = sqrt(dist2);
	cout << "insideLoopAction(), id: " << mol->id() << ", dist2: " << dist2 << ", dist2_min: " << dist2_min << ", dist:" << dist << endl;
	double invDist = 1./dist;
	double fac = invDist * _dReplacement;
	cout << "fac: " << fac << endl;
	for(uint8_t d=0; d<3; ++d)
	{
		double old = mol->r(d);
		double repl = distVec[d] * fac;
		cout << "old: " << old << "repl: " << repl << endl;
		mol->setr(d, old+repl);
	}

	// DEBUG
	dist2 = 0.0;
	for(uint8_t d=0; d<3; ++d)
	{
		distVec[d] = mol->r(d) - _selectedMoleculePos.at(d);
		dist2 += distVec[d]*distVec[d];
	}
	dist = sqrt(dist2);
	cout << "insideLoopAction(), id: " << mol->id() << ", after replacement: " << _dReplacement << ", dist:" << dist << endl;
	// DEBUG

	_numReplacements++;
}

void BubbleMethod::postLoopAction()
{
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();
	int rank = domainDecomp.getRank();
	cout << "rank[" << rank << "]: postLoopAction(), state:" << _nState << endl;
	ParticleContainer* particleCont = global_simulation->getMolecules();

	uint32_t numReplacementsGlobal;
	MPI_Allreduce( &_numReplacements, &numReplacementsGlobal, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
	cout << "rank[" << rank << "]: numReplacementsGlobal=" << numReplacementsGlobal << endl;

	if( !(BMS_GROWING_BUBBLE == _nState && numReplacementsGlobal == 0) )
		return;

	Component* comp2 = global_simulation->getEnsemble()->getComponent(1);  // H2
	Component* comp6 = global_simulation->getEnsemble()->getComponent(5);  // H2 insert

	ParticleContainer* particles = global_simulation->getMolecules();
	ParticleIterator pit;
	for( pit  = particles->iteratorBegin();
			pit != particles->iteratorEnd();
		 ++pit )
	{
		if(pit->id() == _selectedMoleculeID) //componentid()+1 == 3)
		{
			pit->setComponent(comp6);
			pit->setr(0, pit->r(0)-_dBubbleRadius*0.5);
			pit->setv(0, -_vm);
		}
		else if(pit->componentid()+1 == 4)
			pit->setComponent(comp2);
	}

	if(rank == _insRank)
	{
		double rx, ry, rz;
		rx = _selectedMoleculePos.at(0);
		ry = _selectedMoleculePos.at(1);
		rz = _selectedMoleculePos.at(2);
		Molecule mol1(++_maxID, comp6, rx+_dBubbleRadius*0.5, ry, rz, _vm, 0, 0);
		Molecule mol2(++_maxID, comp6, rx, ry+_dBubbleRadius*0.5, rz, 0,  _vm, 0);
		Molecule mol3(++_maxID, comp6, rx, ry-_dBubbleRadius*0.5, rz, 0, -_vm, 0);
		Molecule mol4(++_maxID, comp6, rx, ry, rz+_dBubbleRadius*0.5, 0, 0,  _vm);
		Molecule mol5(++_maxID, comp6, rx, ry, rz-_dBubbleRadius*0.5, 0, 0, -_vm);

		cout << "rank[" << rank << "]: postLoopAction() adding particle at: " << rx << ", " << ry << ", " << rz << endl;
		particleCont->addParticle(mol1, true, true, true);
		particleCont->addParticle(mol2, true, true, true);
		particleCont->addParticle(mol3, true, true, true);
		particleCont->addParticle(mol4, true, true, true);
		particleCont->addParticle(mol5, true, true, true);
	}

	// update maxID
	MPI_Bcast( &_maxID, 1, MPI_UNSIGNED_LONG, _insRank, MPI_COMM_WORLD);

	_nState = BMS_IDLE;
}

void BubbleMethod::postEventNewTimestepAction()
{
	if(_nState != BMS_GROWING_BUBBLE)
	{
//		cout << "postEventNewTimestepAction(), state:" << _nState << endl;
		return;
	}

	ParticleContainer* particles = global_simulation->getMolecules();
	ParticleIterator pit;
	for( pit  = particles->iteratorBegin();
			pit != particles->iteratorEnd();
		 ++pit )
	{
		if(pit->id() == _selectedMoleculeID)
		{
			for(uint8_t d=0; d<3; ++d)
				pit->setr(d, _selectedMoleculePos.at(d) );
		}
	}
}

void BubbleMethod::postUpdateForcesAction()
{
	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();
	int rank = domainDecomp.getRank();
	cout << "rank[" << rank << "]: postUpdateForcesAction(), state:" << _nState << endl;

	if(_nState != BMS_GROWING_BUBBLE)
	{
//		cout << "postEventNewTimestepAction(), state:" << _nState << endl;
		return;
	}

	ParticleContainer* particles = global_simulation->getMolecules();
	ParticleIterator pit;
	for( pit  = particles->iteratorBegin();
			pit != particles->iteratorEnd();
		 ++pit )
	{
		// skip selected molecule, but update its position
	//	cout << "_selectedMoleculeID="<<_selectedMoleculeID<<",mol->id()="<<mol->id()<<endl;
		if(pit->id() == _selectedMoleculeID)
		{
	//		for(uint8_t d=0; d<3; ++d)
	//			_selectedMoleculePos.at(d) = pit->r(d);
	//		cout << "insideLoopAction(): molecule: id=" << pit->id() << ", pos: " << pit->r(0) << ", " << pit->r(1) << ", " << pit->->r(2) << endl;
			continue;
		}
	//	else
	//		global_log->error() << "Selected molecule NOT found!" << endl;


		double dist2_min = (_dBubbleRadius+_dMoleculeOuterRadius)*(_dBubbleRadius+_dMoleculeOuterRadius);
		Component* comp2 = global_simulation->getEnsemble()->getComponent(1);  // H2
		Component* comp5 = global_simulation->getEnsemble()->getComponent(4);  // H2 replace
		if(pit->componentid()+1 == 5)
			pit->setComponent(comp2);
		double pbc[7][3] =
		{
			{0,0,0},
			{1,0,0},
			{0,1,0},
			{0,0,1},
			{-1,0,0},
			{0,-1,0},
			{0,0,-1},
		};
		Domain* domain = global_simulation->getDomain();
		double box[3];
		for(uint8_t d=0; d<3; ++d)
			box[d] = domain->getGlobalLength(d);

		bool bInsideRadius = false;
		double dist2;
		double distVec[3];
		for(uint8_t pbi=0; pbi<7; ++pbi)
		{
			dist2 = 0.0;
			for(uint8_t d=0; d<3; ++d)
			{
				distVec[d] = (pit->r(d)+pbc[pbi][d]*box[d]) - _selectedMoleculePos.at(d);
				dist2 += distVec[d]*distVec[d];
			}
	//		cout << "postUpdateForcesAction(), id: " << pit->id() << ", dist2: " << dist2 << ", dist2_min: " << dist2_min << endl;
			if(dist2 < dist2_min)
			{
				bInsideRadius = true;
				break;
			}
		}

		if(false == bInsideRadius)
			continue;
		else
			pit->setComponent(comp5);

		double dist = sqrt(dist2);
		cout << "rank[" << rank << "]: postUpdateForcesAction(), id: " << pit->id() << ", dist2: " << dist2 << ", dist2_min: " << dist2_min << ", dist:" << dist << endl;
		double invDist = 1./dist;
		double Fadd[3];
		cout << "rank[" << rank << "]: postUpdateForcesAction(), id: " << pit->id() << ", F: " << pit->F(0) << ", " << pit->F(1) << ", " << pit->F(2) << endl;
		for(uint8_t d=0; d<3; ++d)
			Fadd[d] = distVec[d] * _dMaxForce * (1-dist2/dist2_min);
		pit->Fadd(Fadd);
		cout << "rank[" << rank << "]: postUpdateForcesAction(), id: " << pit->id() << ", F: " << pit->F(0) << ", " << pit->F(1) << ", " << pit->F(2) << endl;

//		// DEBUG
//		dist2 = 0.0;
//		for(uint8_t d=0; d<3; ++d)
//		{
//			distVec[d] = pit->r(d) - _selectedMoleculePos.at(d);
//			dist2 += distVec[d]*distVec[d];
//		}
//		dist = sqrt(dist2);
//		cout << "insideLoopAction(), id: " << pit->id() << ", after replacement: " << _dReplacement << ", dist:" << dist << endl;
//		// DEBUG

		_numReplacements++;
	}
}
