/*
 * ParticleTracker.h
 *
 *  Created on: 27.01.2017
 *      Author: mheinen
 */

#ifndef PARTICLETRACKER_H_
#define PARTICLETRACKER_H_

#include <vector>
#include <map>
#include <string>
#include "molecules/MoleculeForwardDeclaration.h"

class DomainDecompBase;

class TrackerState;
class ParticleTracker
{

public:
	ParticleTracker(DomainDecompBase* domainDecomp, unsigned long simstepStart,
			unsigned long simstepStop);
	~ParticleTracker();

	void Prepare();
	void PreLoopAction(unsigned long simstep);
	void LoopAction(Molecule* mol);
	void PostLoopAction();

private:
	unsigned int ReadMoleculeIDs();
	void ShowSampleMap();
	void ResetSampleMap();
	void CreateDataFiles();
	void CreateMoleculeIDsFile();
	void InitDatastructures();

	// methods used by state classes
	unsigned long GetSimstep() {return _simstep;}
	void SetSimstep(unsigned long simstep) {_simstep = simstep;}
	void CalcGlobalValuesData();
	void WriteToFileData();
	void WriteToFileIDs();
	void ClearSampleDataVector() {_vecSampleDataDouble.clear();}
	unsigned long GetWriteFreqData() {return _writeFreqData;}
	void AppendDataMap() {_vecSampleDataDouble.push_back(_mapMolIDsDoubleVals);}
	void CollectIDs(Molecule* mol);
	void SampleData(Molecule* mol);
	bool SimstepInsideRange() {return _simstep >= _simstepStart && _simstep < _simstepStop;}

	void ChangeState(TrackerState* state) {_trackerState = state;}

private:

	friend class TrackerState;
	friend class TStateSleep;
	friend class TStateCollectMolIDs;
	friend class TStateCollectData;

	DomainDecompBase* _domainDecomp;
	std::vector<unsigned long> _vecMoleculeIDs;
	std::vector<std::string> _vecQuantyNamesDouble;
	std::vector<std::map<unsigned long, std::vector<double> > > _vecSampleDataDouble;
	std::map<unsigned long, std::vector<double> > _mapMolIDsDoubleVals;
	std::vector<std::map<unsigned long, unsigned char> > _mapCompIDs;
	std::map<unsigned long, char> _mapCollectIDs;

	unsigned long _simstep;

	TrackerState* _trackerState;
	unsigned long _writeFreqData;
	unsigned long _simstepStart;
	unsigned long _simstepStop;
};


class ParticleTracker;
class TrackerState
{
protected:
	TrackerState() {}

public:
	virtual void PreLoopAction(ParticleTracker* context, unsigned long simstep) = 0;
	virtual void LoopAction(ParticleTracker* context, Molecule* mol) = 0;
	virtual void PostLoopAction(ParticleTracker* context) = 0;

protected:
	void ChangeState(ParticleTracker* tracker, TrackerState* state);
};

class TStateSleep : public TrackerState
{
public:
	static TrackerState* Exemplar();

protected:
	TStateSleep() {}

	virtual void PreLoopAction(ParticleTracker* context, unsigned long simstep);
	virtual void LoopAction(ParticleTracker* context, Molecule* mol);
	virtual void PostLoopAction(ParticleTracker* context);

private:
	static TrackerState* _exemplar;
};

class TStateCollectMolIDs : public TrackerState
{
public:
	static TrackerState* Exemplar();

protected:
	TStateCollectMolIDs() {}

	virtual void PreLoopAction(ParticleTracker* context, unsigned long simstep);
	virtual void LoopAction(ParticleTracker* context, Molecule* mol);
	virtual void PostLoopAction(ParticleTracker* context);

private:
	static TrackerState* _exemplar;
};

class TStateCollectData : public TrackerState
{
public:
	static TrackerState* Exemplar();

protected:
	TStateCollectData() {}

	virtual void PreLoopAction(ParticleTracker* context, unsigned long simstep);
	virtual void LoopAction(ParticleTracker* context, Molecule* mol);
	virtual void PostLoopAction(ParticleTracker* context);

private:
	static TrackerState* _exemplar;
};


#endif /* PARTICLETRACKER_H_ */
