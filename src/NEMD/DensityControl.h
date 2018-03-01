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
#include "utils/CommVar.h"
#include "utils/CtrlVar.h"

#include <vector>
#include <list>
#include <unordered_map>
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
class ParticleManipDirector;
class MainLoopAction;

namespace dec
{

struct CompVarsStruct
{
	CtrlVar<CommVar<double> > density;
	CtrlVar<CommVar<int64_t> > numMolecules;
	uint32_t compID;
	uint32_t proxyID;
	std::list<uint64_t> particleIDs;
};

enum ControlRegionStates : int32_t
{
	CRS_UNKNOWN = -1,
	CRS_BALANCED = 0,
	CRS_CHANGE_MOLECULES = 1,
	CRS_INSERT_MOLECULES = 2,
};

typedef std::unordered_map<uint32_t, CommVar<uint64_t> > cid_count_map;
struct manip_count_struct
{
	CommVar<uint64_t> deleted;
	CommVar<uint64_t> inserted;
	cid_count_map changed_to_from;
};
typedef std::unordered_map<uint32_t, manip_count_struct> cid_manip_count_map;
typedef std::unordered_map<uint32_t, manip_count_struct>::iterator cid_manip_count_map_it;

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

    void MeasureDensity(Simulation* simulation, Molecule* mol);
    void CalcGlobalValues(Simulation* simulation);
    void ManipulateParticles(Simulation* simulation, Molecule* mol);
    void ManipulateParticleForces(Simulation* simulation, Molecule* mol);
    void FinalizeParticleManipulation(Simulation* simulation, MainLoopAction* action);

    void ResetLocalValues(Simulation* simulation);
    void UpdateVolume(DomainDecompBase* domainDecomp);

    void WriteHeaderDeletedMolecules();
    void WriteDataDeletedMolecules(unsigned long simstep);
	void WriteHeaderParticleManipCount();
	void WriteDataParticleManipCount(unsigned long simstep);

	// Connection to MettDeamon
	void ConnectMettDeamon(const std::vector<MettDeamon*>& mettDeamon) {_mettDeamon = _bMettDeamonConnected ? mettDeamon.at(_nMettDeamonInstanceIndex) : NULL;}

	// getters
	uint32_t getGlobalMinNumMoleculesSpreadCompID();
	uint32_t getGlobalMaxNumMoleculesSpreadCompID();
	void getGlobalMinMaxNumMoleculesSpreadCompIDs(uint32_t& cidMinSpread, uint32_t& cidMaxSpread);
	std::vector<CompVarsStruct> getCompVars() {return _compVars;}
	std::list<uint64_t> GetLocalParticleIDs(const uint32_t& nCompID) {return _compVars.at(nCompID).particleIDs;}
	int64_t getLocalNumMoleculesSpread(uint32_t nCompID) {return _compVars.at(nCompID).numMolecules.spread.local;}
	bool getVacuum() {return _bVacuum;}

	// checks
	bool globalTargetDensityExeeded(uint32_t cid)   {return _compVars.at(cid).numMolecules.spread.global > 0;}
	bool globalTargetDensityUndershot(uint32_t cid) {return _compVars.at(cid).numMolecules.spread.global < 0;}
	bool globalCompositionBalanced();

	// inform about particle manipulations
	void informParticleInserted(Molecule mol);
	void informParticleDeleted(Molecule mol);
	void informParticleChanged(Molecule from, Molecule to);

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

private:
	bool _bVacuum;

	// parameter
	std::vector<CompVarsStruct> _compVars;

	CommVar<double> _dVolume;
	CommVar<double> _dInvertVolume;

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
	ParticleManipDirector* _director;

	cid_manip_count_map _manip_count_map;
};

}  // namespace dec

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
    void ResetLocalValues(Simulation* simulation);
    void MeasureDensity(Simulation* simulation, Molecule* mol);
    void CalcGlobalValues(Simulation* simulation);
    void ManipulateParticles(Simulation* simulation, Molecule* mol);
    void ManipulateParticleForces(Simulation* simulation, Molecule* mol);
    void FinalizeParticleManipulation(Simulation* simulation, MainLoopAction* action);

    unsigned long GetControlFreq() {return _nControlFreq;}
    unsigned long GetStart() {return _nStart;}
    unsigned long GetStop()  {return _nStop;}

    bool ProcessIsRelevant() {return _bProcessIsRelevant;}

    void WriteHeaderDeletedMolecules();
    void WriteDataDeletedMolecules(Simulation* simulation);
	void WriteHeaderParticleManipCount();
	void WriteDataParticleManipCount(Simulation* simulation);

	// NEMD flags
	void     SetFlagsNEMD(const uint32_t &flagsNEMD) {_flagsNEMD = _flagsNEMD | flagsNEMD;}
	uint32_t GetFlagsNEMD() {return _flagsNEMD;}
	bool CheckFlagNEMD(uint32_t nFlag) {return (_flagsNEMD & nFlag);}

	// Connection to MettDeamon
	void ConnectMettDeamon(const std::vector<MettDeamon*>& mettDeamon);

	MainLoopAction* getPreForceAction() {return _preForceAction;}
	MainLoopAction* getPostForceAction() {return _postForceAction;}

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
	virtual void preFirstLoop(Simulation* simulation) = 0;
	virtual void insideFirstLoop(Simulation* simulation, Molecule* mol) = 0;
	virtual void postFirstPreSecondLoop(Simulation* simulation) = 0;
	virtual void insideSecondLoop(Simulation* simulation, Molecule* mol) = 0;
	virtual void postSecondLoop(Simulation* simulation) = 0;
};

class PreForceAction : public MainLoopAction
{
public:
	PreForceAction(DensityControl* parent) : MainLoopAction(parent) {}
	virtual ~PreForceAction() {}

	virtual void preFirstLoop(Simulation* simulation);
	virtual void insideFirstLoop(Simulation* simulation, Molecule* mol);
	virtual void postFirstPreSecondLoop(Simulation* simulation);
	virtual void insideSecondLoop(Simulation* simulation, Molecule* mol);
	virtual void postSecondLoop(Simulation* simulation);
};

class PostForceAction : public MainLoopAction
{
public:
	PostForceAction(DensityControl* parent) : MainLoopAction(parent) {}
	virtual ~PostForceAction() {}

	virtual void preFirstLoop(Simulation* simulation);
	virtual void insideFirstLoop(Simulation* simulation, Molecule* mol);
	virtual void postFirstPreSecondLoop(Simulation* simulation);
	virtual void insideSecondLoop(Simulation* simulation, Molecule* mol);
	virtual void postSecondLoop(Simulation* simulation);
};

#endif /* DENSITYCONTROL_H_ */
