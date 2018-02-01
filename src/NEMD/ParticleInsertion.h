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
#include "utils/CommVar.h"
#include "utils/CtrlVar.h"

template<typename T1, typename T2>
void select_rnd_elements(std::list<T1>& mylist, std::vector<T1>& myvec, T2 numSelect);

template<typename T>
class ChangeVar
{
public:
	T from;
	T to;
};

namespace dec
{
	class ControlRegion;
	struct CompVarsStruct;
}
class MainLoopAction;
class Simulation;
class XMLfileUnits;
class ParticleManipulator;
class ParticleDeleter;
class ParticleInsertion;
class ParticleManipDirector
{
public:
	ParticleManipDirector(dec::ControlRegion* region);
	~ParticleManipDirector() {}

	void readXML(XMLfileUnits& xmlconfig);
	void localValuesReseted(Simulation* simulation);
	void globalValuesCalculated(Simulation* simulation);
	void ManipulateParticles(Simulation* simulation, Molecule* mol);
	void ManipulateParticleForces(Simulation* simulation, Molecule* mol);
	void FinalizeParticleManipulation(Simulation* simulation, MainLoopAction* action);
	std::vector<dec::CompVarsStruct> getCompVars();
	double GetLowerCorner(uint32_t nDim);
	double GetWidth(uint32_t nDim);
	ChangeVar<uint32_t> getNextChangeIDs() {return _nextChangeIDs;}
	std::list<uint64_t> GetLocalParticleIDs(const uint32_t& nCompID);
	int64_t getLocalNumMoleculesSpread(uint32_t nCompID);

private:
	bool setNextChangeIDs();
	bool setNextInsertIDs();

private:
	dec::ControlRegion* _region;
	ParticleDeleter* _deleter;
	ParticleManipulator* _manipulator;
	std::vector<ParticleInsertion*> _inserter;
	std::vector<ParticleInsertion*> _changer;
	ChangeVar<uint32_t> _nextChangeIDs;
	MainLoopAction* _preForceAction;
	MainLoopAction* _postForceAction;
};

enum ParticleManipulatorStates : uint32_t
{
	PMS_IDLE = 1,
	PMS_BUSY = 2,
};

class ParticleManipulator
{
protected:
	ParticleManipulator(ParticleManipDirector* director) : _director(director), _nManipState(PMS_IDLE) {}
	virtual ~ParticleManipulator() {}

public:
	virtual void readXML(XMLfileUnits& xmlconfig) = 0;
	virtual void Reset(Simulation* simulation) = 0;
	virtual void PrepareParticleManipulation(Simulation* simulation) = 0;
	virtual void ManipulateParticles(Simulation* simulation, Molecule* mol) = 0;
	virtual void ManipulateParticleForces(Simulation* simulation, Molecule* mol) = 0;
	virtual void FinalizeParticleManipulation(Simulation* simulation) = 0;
	virtual void FinalizeParticleManipulation_preForce(Simulation* simulation) = 0;

	CommVar<uint64_t> getNumManipulatedParticles() {return _numManipulatedParticles;}
	uint32_t getManipState() {return _nManipState;}

protected:
	ParticleManipDirector* _director;
	CommVar<uint64_t> _numManipulatedParticles;
	uint32_t _nManipState;
};

class ParticleDeleter : public ParticleManipulator
{
public:
	ParticleDeleter(ParticleManipDirector* director);
	virtual ~ParticleDeleter() {}

	virtual void readXML(XMLfileUnits& xmlconfig) {}
	virtual void Reset(Simulation* simulation) {}
	virtual void PrepareParticleManipulation(Simulation* simulation);
	virtual void ManipulateParticles(Simulation* simulation, Molecule* mol);
	virtual void ManipulateParticleForces(Simulation* simulation, Molecule* mol) {}
	virtual void FinalizeParticleManipulation(Simulation* simulation) {_nManipState = PMS_IDLE;}
	virtual void FinalizeParticleManipulation_preForce(Simulation* simulation) {_nManipState = PMS_IDLE;}
	void CreateDeletionLists(std::vector<dec::CompVarsStruct> compVars);

private:
	std::vector< std::vector<uint64_t> > _deletionLists;
};

class Simulation;
class CuboidRegion;
class ParticleInsertion : public ParticleManipulator
{
protected:
	ParticleInsertion(ParticleManipDirector* director, uint32_t state);
	~ParticleInsertion() {}

public:
	virtual void Reset(Simulation* simulation) = 0;
	virtual void PrepareParticleManipulation(Simulation* simulation) = 0;
	virtual void ManipulateParticles(Simulation* simulation, Molecule* mol) = 0;
	virtual void ManipulateParticleForces(Simulation* simulation, Molecule* mol) = 0;
	virtual void FinalizeParticleManipulation(Simulation* simulation) = 0;
	virtual void FinalizeParticleManipulation_preForce(Simulation* simulation) = 0;

	// getters / setters
	void setState(uint32_t state) {_nState = state;}
	uint32_t getState() {return _nState;}
	uint32_t getTargetCompID() {return _nTargetCompID;}
	void setTargetCompID(uint32_t nVal) {_nTargetCompID = nVal;}

protected:
	uint32_t _nState;
	uint32_t _nTargetCompID;
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

struct BubbleVars
{
	CtrlVar<CommVar<double> > radius;
	CtrlVar<CommVar<double> > radiusSquared;
	double velocity;

	struct ForceVars
	{
		double radius;
		double maxVal;
	} force;
};

struct ID_pair
{
	uint64_t mid;
	uint32_t cid;
};

enum BubbleMethodTypes
{
	BMT_CHANGER = 1,
	BMT_INSERTER = 2,
};

class ParticleSelector;
class XMLfileUnits;
class Random;
class BubbleMethod : public ParticleInsertion
{
public:
	BubbleMethod(ParticleManipDirector* director, uint32_t nType);
	~BubbleMethod();

	virtual void readXML(XMLfileUnits& xmlconfig);
	virtual void Reset(Simulation* simulation);
	virtual void PrepareParticleManipulation(Simulation* simulation);
	virtual void ManipulateParticles(Simulation* simulation, Molecule* mol);
	virtual void ManipulateParticleForces(Simulation* simulation, Molecule* mol);
	virtual void FinalizeParticleManipulation(Simulation* simulation);
	virtual void FinalizeParticleManipulation_preForce(Simulation* simulation);

	// Selector info
	double GetLowerCorner(uint32_t nDim);
	double GetWidth(uint32_t nDim);
	std::list<uint64_t> GetLocalParticleIDs(const uint32_t& nCompID);
	int64_t getLocalNumMoleculesSpread(uint32_t nCompID);

	// getters / setters
	uint32_t getType() {return _nType;}
	ChangeVar<uint32_t> getNextChangeIDs() {return _nextChangeIDs;}

private:
	void selectParticle(Simulation* simulation);
	void initInsertionMolecules(Simulation* simulation);  // inserter only
	void FreezeSelectedMolecule(Molecule* mol);
	void resetBubbleMoleculeComponent(Molecule* mol);
	double calcMinSquaredDistance(Molecule* mol);
	bool outerMoleculeRadiusCutsBubbleRadius(Simulation* simulation, Molecule* mol, const double& dist2_min, double& dist2, double* distVec);
	void updateActualBubbleRadiusLocal(Molecule* mol, const double& dist);
	void updateForceOnBubbleCuttingMolecule(Molecule* mol, const double& dist2_min, const double& dist2, const double& dist, double* distVec);
	void GrowBubble(Simulation* simulation, Molecule* mol);
	void ChangeIdentity(Simulation* simulation, Molecule* mol);
	bool checkBubbleRadius();

private:
	ParticleSelector* _selector;
	BubbleVars _bubble;
	int _insRank;
	uint64_t _maxID;
	uint64_t _selectedMoleculeID;
	uint32_t _nType;
	std::array<double, 3> _daInsertionPosition;
	std::array<double, 3> _selectedMoleculeInitPos;
	ChangeVar<uint32_t> _nextChangeIDs;
	struct {
		std::vector<Molecule> initial;
		std::vector<Molecule> actual;
	} _insertMolecules;
	std::vector<ID_pair> _bubbleMolecules;
	struct {
		std::vector<std::array<double, 3> > initial;
		std::vector<std::array<double, 3> > actual;
	} _pbc;
};

class ParticleSelector
{
protected:
	ParticleSelector(BubbleMethod* parent) : _parent(parent){}
	virtual ~ParticleSelector() {}

public:
	virtual int selectParticle(Simulation* simulation, Molecule& mol) = 0;
protected:
	virtual void collectInfo() = 0;
protected:
	BubbleMethod* _parent;
};

class RandomSelector : public ParticleSelector
{
public:
	RandomSelector(BubbleMethod* parent);
	virtual ~RandomSelector() {}

public:
	virtual int selectParticle(Simulation* simulation, Molecule& mol);
protected:
	virtual void collectInfo();
private:
	void generateRandomInsPos();
private:
	Random* _random;
	std::array<double, 3> _daInsertionPosition;
	std::array<double, 3> _insRegLowerCorner;
	std::array<double, 3> _insRegWidth;
};

class CompDependSelector : public ParticleSelector
{
public:
	CompDependSelector(BubbleMethod* parent);
	virtual ~CompDependSelector() {}

public:
	virtual int selectParticle(Simulation* simulation, Molecule& mol);
protected:
	virtual void collectInfo();
};

#endif /* PARTICLEINSERTION_H_ */
