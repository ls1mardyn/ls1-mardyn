/*
 * ParticleInsertion.h
 *
 *  Created on: 12.01.2018
 *      Author: mheinen
 */

#ifndef PARTICLEINSERTION_H_
#define PARTICLEINSERTION_H_

#include <cstdint>
#include <array>
#include <vector>

#include "molecules/Molecule.h"

class CuboidRegion;
class ParticleInsertion
{
protected:
	ParticleInsertion(CuboidRegion* parent, uint32_t state);
	~ParticleInsertion() {}

public:
	virtual void preLoopAction() = 0;
	virtual void insideLoopAction(Molecule* mol) = 0;
	virtual void postLoopAction() = 0;

	void setState(uint32_t state) {_nState = state;}

protected:
	CuboidRegion* _parent;
	uint32_t _nState;
};

enum BubbleMethodStates
{
	BMS_IDLE = 0,
	BMS_SELECT_RANDOM_SPOT = 1,
	BMS_GROWING_BUBBLE = 2,
	BMS_REACHED_BUBBLE_SIZE = 3,
};

class XMLfileUnits;
class Random;
class BubbleMethod : public ParticleInsertion
{
public:
	BubbleMethod(CuboidRegion* parent);
	~BubbleMethod();

	virtual void preLoopAction();
	virtual void insideLoopAction(Molecule* mol);
	virtual void postLoopAction();

	void readXML(XMLfileUnits& xmlconfig);

private:
	Random* _random;
	double _dBubbleForce;
	double _dBubbleRadius;
	double _dBubbleRadiusSquared;
	double _dMoleculeOuterRadius;
	double _dReplacement;
	uint32_t _numReplacements;
	uint64_t _maxID;
	std::array<double, 3> _daInsertionPosition;

};

#endif /* PARTICLEINSERTION_H_ */
