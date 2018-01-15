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
	_numReplacements(0),
	_maxID(100000000)
{
	_random = new Random();

	// tmp
	_dMoleculeOuterRadius = 4.3675 * 0.5;
	_dBubbleRadius = _dMoleculeOuterRadius;
	_dReplacement = _dMoleculeOuterRadius * 2 * 0.01;
}

BubbleMethod::~BubbleMethod()
{
	delete _random;
}

void BubbleMethod::preLoopAction()
{
	cout << "preLoopAction(), state:" << _nState << endl;
	// Reset replacement counter
	_numReplacements = 0;

	if(_nState != BMS_SELECT_RANDOM_SPOT)
		return;

	DomainDecompBase domainDecomp = global_simulation->domainDecomposition();
	Domain* domain = global_simulation->getDomain();
	double bbMin[3];
	double bbMax[3];
	domainDecomp.getBoundingBoxMinMax(domain, bbMin, bbMax);

	for(uint8_t d=0; d<3; ++d)
	{
		float rnd = _random->rnd();
		cout << "rnd=" << rnd << endl;
//		_daInsertionPosition.at(d) = rnd * (bbMax[d] - bbMin[d]);
		_daInsertionPosition.at(d) = _parent->GetLowerCorner(d) + rnd * _parent->GetWidth(d);
	}
	// DEBUG
	_daInsertionPosition.at(2) = global_simulation->getDomain()->getGlobalLength(2) * 0.5;
	// DEBUG

	_nState = BMS_GROWING_BUBBLE;
}

void BubbleMethod::insideLoopAction(Molecule* mol)
{
	if(_nState != BMS_GROWING_BUBBLE)
		return;

	double dist2 = 0.0;
	double distVec[3];
	for(uint8_t d=0; d<3; ++d)
	{
		distVec[d] = _daInsertionPosition.at(d) - mol->r(d);
		dist2 += distVec[d]*distVec[d];
	}

	if(dist2 > (_dBubbleRadius+_dMoleculeOuterRadius)*(_dBubbleRadius+_dMoleculeOuterRadius) )
		return;

	double dist = sqrt(dist2);
	cout << "insideLoopAction(), dist:" << dist << endl;
	double invDist = 1./dist;
	for(uint8_t d=0; d<3; ++d)
	{
		double repl = distVec[d] * invDist * _dReplacement;
		mol->setr(d, repl);
	}
	_numReplacements++;
}

void BubbleMethod::postLoopAction()
{
	cout << "postLoopAction(), state:" << _nState << endl;
	ParticleContainer* particleCont = global_simulation->getMolecules();

	if(BMS_GROWING_BUBBLE == _nState && _numReplacements > 0)
		return;

	double rx, ry, rz;
	rx = _daInsertionPosition.at(0);
	ry = _daInsertionPosition.at(1);
	rz = _daInsertionPosition.at(2);
	Component* comp3 = global_simulation->getEnsemble()->getComponent(2);
	Molecule mol(++_maxID, comp3, rx, ry, rz);

	cout << "postLoopAction() adding particle at: " << rx << ", " << ry << ", " << rz << endl;
	particleCont->addParticle(mol, true, false, true);

	_nState = BMS_IDLE;
}

void BubbleMethod::readXML(XMLfileUnits& xmlconfig)
{

}


