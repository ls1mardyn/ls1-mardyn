/*
 * DensityControl.h
 *
 *  Created on: 29.05.2015
 *      Author: mheinen
 */

#ifndef DENSITYCONTROL_H_
#define DENSITYCONTROL_H_

#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "utils/ObserverBase.h"
#include "utils/Region.h"
#include "Common.h"

#include <vector>
#include <list>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

class MettDeamon;
class XMLfileUnits;
class Simulation;
class Domain;
class DomainDecompBase;
class DensityControl;
class ParticleInsertion;

namespace dec
{

struct compVars
{
	controlVar<double> density;
	controlVar<int64_t> numMolecules;
	uint32_t compID;
	uint32_t proxyID;
	ParticleInsertion* insertion;
	std::list<uint64_t> particleIDs;
	std::vector<uint64_t> deletionList;
};

enum ControlRegionStates : int32_t
{
	CRS_UNKNOWN = -1,
	CRS_BALANCED = 0,
	CRS_CHANGE_MOLECULES = 1,
	CRS_INSERT_MOLECULES = 2,
};

class ControlRegion : public CuboidRegionObs
{
public:
	ControlRegion(DensityControl* parent, double dLowerCorner[3], double dUpperCorner[3] );
	ControlRegion(DensityControl* parent, double dLowerCorner[3], double dUpperCorner[3], unsigned int nTargetComponentID, const double& dTargetDensity);
	virtual ~ControlRegion();

	void readXML(XMLfileUnits& xmlconfig);
    void CheckBounds();
    void Init(DomainDecompBase* domainDecomp);
    void InitMPI();
    bool ProcessIsRelevant() {return _bProcessIsRelevant;}

    double GetDensityGlobal(uint8_t cid) {return _compVars.at(cid).density.actual.global;}

    void MeasureDensity(Molecule* mol);
    void CalcGlobalValues();
    void CreateDeletionLists();

    void ControlDensity(Molecule* mol, Simulation* simulation);
    void postLoopAction();
    void postEventNewTimestepAction();
    void postUpdateForcesAction();

    void ResetLocalValues();
    void UpdateVolume(DomainDecompBase* domainDecomp);

    void WriteHeaderDeletedMolecules();
    void WriteDataDeletedMolecules(unsigned long simstep);

	// Connection to MettDeamon
	void ConnectMettDeamon(const std::vector<MettDeamon*>& mettDeamon) {_mettDeamon = _bMettDeamonConnected ? mettDeamon.at(_nMettDeamonInstanceIndex) : NULL;}
private:
	void CheckState()
	{
		if(_compVars.at(0).numMolecules.actual.global < _compVars.at(0).numMolecules.target.global)
			_nState = CRS_INSERT_MOLECULES;
		else
			_nState = CRS_BALANCED;
		for(auto comp : _compVars)
		{
			if(comp.numMolecules.actual.global < comp.numMolecules.target.global)
			{
				_nState = CRS_CHANGE_MOLECULES;
				break;
			}
		}
	}
	uint32_t NextInsertionID()
	{
		uint32_t cid = 0;
		uint64_t max = 0;
		for(auto comp : _compVars)
		{
			if(0 == comp.compID)
				continue;

			uint64_t diff = comp.numMolecules.target.global - comp.numMolecules.actual.global;
			if(diff > max)
			{
				max = diff;
				cid = comp.compID;
			}
		}
		return cid;
	}
	bool InsertionAllIdle();

	template<typename T1, typename T2>
	void select_rnd_elements(std::list<T1>& mylist, std::vector<T1>& myvec, T2 numSelect);

private:
	// parameter
	std::vector<compVars> _compVars;

	commVar<double> _dVolume;
	commVar<double> _dInvertVolume;

    int* _ranks;
    bool _bProcessIsRelevant;
    MPI_Comm _newcomm;

    // deleted molecules data
    unsigned long _nDeletedNumMoleculesLocal;
    unsigned long _nDeletedNumMoleculesGlobal;
    double _dDeletedEkinLocal[3];
    double _dDeletedEkinGlobal[3];
    double _dDeletedVelocityLocal[3];
    double _dDeletedVelocityGlobal[3];
    double _dDeletedVelocitySquaredLocal[3];
    double _dDeletedVelocitySquaredGlobal[3];

    // instances / ID
    static unsigned short _nStaticID;

	// Connection to MettDeamon
	MettDeamon* _mettDeamon;
	bool _bMettDeamonConnected;
	uint8_t _nMettDeamonInstanceIndex;

	// identity change (by component ID)
	std::vector<uint32_t> _vecChangeCompIDs;

	int32_t _nState;
};

}

class MainLoopAction;
class DensityControl : public ControlInstance
{
public:
	DensityControl(DomainDecompBase* domainDecomp, Domain* domain);
    DensityControl(DomainDecompBase* domainDecomp, Domain* domain, unsigned long nControlFreq, unsigned long nStart, unsigned long nStop);
	virtual ~DensityControl();

	void readXML(XMLfileUnits& xmlconfig);
    std::string GetShortName() {return "DeC";}
    void AddRegion(dec::ControlRegion* region);
    int GetNumRegions() {return _vecControlRegions.size();}
    dec::ControlRegion* GetControlRegion(unsigned short nRegionID) {return _vecControlRegions.at(nRegionID-1); }  // vector index starts with 0, region index with 1

	void preForce_action(Simulation* simulation);
	void postForce_action(Simulation* simulation);

    void CheckRegionBounds();
    void Init(DomainDecompBase* domainDecomp);
    void InitMPI();
    void ResetLocalValues(unsigned long simstep);
    void MeasureDensity(Molecule* mol, unsigned long simstep);
    void CalcGlobalValues(unsigned long simstep);
    void CreateDeletionLists();
    void ControlDensity(Molecule* mol, Simulation* simulation);
    void postLoopAction();
    void postEventNewTimestepAction();
    void postUpdateForcesAction();

    unsigned long GetControlFreq() {return _nControlFreq;}
    unsigned long GetStart() {return _nStart;}
    unsigned long GetStop()  {return _nStop;}

    bool ProcessIsRelevant() {return _bProcessIsRelevant;}

    void WriteHeaderDeletedMolecules();
    void WriteDataDeletedMolecules(unsigned long simstep);

	// NEMD flags
	void     SetFlagsNEMD(const uint32_t &flagsNEMD) {_flagsNEMD = _flagsNEMD | flagsNEMD;}
	uint32_t GetFlagsNEMD() {return _flagsNEMD;}
	bool CheckFlagNEMD(uint32_t nFlag) {return (_flagsNEMD & nFlag);}

	// Connection to MettDeamon
	void ConnectMettDeamon(const std::vector<MettDeamon*>& mettDeamon);

private:
    std::vector<dec::ControlRegion*> _vecControlRegions;
	uint64_t _nStart;
	uint64_t _nControlFreq;
	uint64_t _nStop;
	uint64_t _nWriteFreqDeleted;

    bool _bProcessIsRelevant;

	uint32_t _flagsNEMD;

	// template methods
	MainLoopAction* _preForceAction;
	MainLoopAction* _postForceAction;
};

class MainLoopAction
{
protected:
	MainLoopAction(DensityControl* parent) : _parent(parent) {}
	virtual ~MainLoopAction() {}

protected:
	DensityControl* _parent;

public:
	void performAction(Simulation* simulation);
	virtual void preFirstLoop(unsigned long simstep) = 0;
	virtual void insideFirstLoop(Molecule* mol, unsigned long simstep) = 0;
	virtual void postFirstPreSecondLoop(unsigned long simstep) = 0;
	virtual void insideSecondLoop(Molecule* mol, Simulation* simulation) = 0;
	virtual void postSecondLoop(unsigned long simstep) = 0;
};

class PreForceAction : public MainLoopAction
{
public:
	PreForceAction(DensityControl* parent) : MainLoopAction(parent) {}
	virtual ~PreForceAction() {}

	virtual void preFirstLoop(unsigned long simstep);
	virtual void insideFirstLoop(Molecule* mol, unsigned long simstep);
	virtual void postFirstPreSecondLoop(unsigned long simstep);
	virtual void insideSecondLoop(Molecule* mol, Simulation* simulation);
	virtual void postSecondLoop(unsigned long simstep);
};

class PostForceAction : public MainLoopAction
{
public:
	PostForceAction(DensityControl* parent) : MainLoopAction(parent) {}
	virtual ~PostForceAction() {}

	virtual void preFirstLoop(unsigned long simstep);
	virtual void insideFirstLoop(Molecule* mol, unsigned long simstep);
	virtual void postFirstPreSecondLoop(unsigned long simstep);
	virtual void insideSecondLoop(Molecule* mol, Simulation* simulation);
	virtual void postSecondLoop(unsigned long simstep);
};

#endif /* DENSITYCONTROL_H_ */
