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

	virtual void readXML(XMLfileUnits& xmlconfig) = 0;

	virtual void preLoopAction() = 0;
	virtual void insideLoopAction(Molecule* mol) = 0;
	virtual void postLoopAction() = 0;
	virtual void postEventNewTimestepAction() = 0;
	virtual void postUpdateForcesAction() = 0;

	void setState(uint32_t state) {_nState = state;}
	uint32_t getState() {return _nState;}

protected:
	CuboidRegion* _parent;
	uint32_t _nState;
};

enum BubbleMethodStates
{
	BMS_IDLE = 0,
	BMS_SELECT_RANDOM_MOLECULE = 1,
	BMS_GROWING_BUBBLE = 2,
	BMS_REACHED_BUBBLE_SIZE = 3,
};

enum InsertionComponentTypes
{
	ICID_ORIGIN = 0,
	ICID_SELECTED = 1,
	ICID_REPLACE = 2,
	ICID_INSERTED = 3,
};

struct InsertionComponentIDs
{
	uint32_t origin;
	uint32_t selected;
	uint32_t replace;
	uint32_t inserted;
};

class XMLfileUnits;
class Random;
class BubbleMethod : public ParticleInsertion
{
public:
	BubbleMethod(CuboidRegion* parent);
	~BubbleMethod();

	virtual void readXML(XMLfileUnits& xmlconfig);

	virtual void preLoopAction();
	virtual void insideLoopAction(Molecule* mol);
	virtual void postLoopAction();
	virtual void postEventNewTimestepAction();
	virtual void postUpdateForcesAction();

private:
	Random* _random;
	double _dBubbleForce;
	double _dBubbleRadius;
	double _dBubbleRadiusSquared;
	double _dMoleculeOuterRadius;
	double _dReplacement;
	double _dMaxForce;
	double _vm;
	int _insRank;
	uint32_t _numReplacements;
	uint64_t _maxID;
	uint64_t _selectedMoleculeID;
	std::array<double, 3> _daInsertionPosition;
	std::array<double, 3> _selectedMoleculePos;

};

#endif /* PARTICLEINSERTION_H_ */
