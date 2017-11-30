/*
 * MettDeamon.h
 *
 *  Created on: 03.04.2017
 *      Author: thet
 */

#ifndef METTDEAMON_H_
#define METTDEAMON_H_

#include "molecules/Molecule.h"
#include "Domain.h"

#include <map>
#include <array>
#include <fstream>
#include <vector>
#include <cstdint>
#include <limits>

#define FORMAT_SCI_MAX_DIGITS std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)

enum ReadReservoirMethods : uint8_t
{
	RRM_UNKNOWN = 0,
	RRM_READ_FROM_FILE = 1,
	RRM_READ_FROM_MEMORY = 2,
	RRM_AMBIGUOUS = 3,
};

enum MovingDirections : uint8_t
{
	MD_UNKNOWN = 0,
	MD_LEFT_TO_RIGHT = 1,
	MD_RIGHT_TO_LEFT = 2,
};

enum FeedRateMethod : uint8_t
{
	FRM_UNKNOWN = 0,
	FRM_DELETED_MOLECULES = 1,
	FRM_CHANGED_MOLECULES = 2,
	FRM_DENSITY = 3
};

class Domain;
class Ensemble;
class DomainDecompBase;
class ParticleContainer;
class XMLfileUnits;

class MettDeamon
{
public:
	MettDeamon();
	~MettDeamon();

	void readXML(XMLfileUnits& xmlconfig);

	uint64_t getnNumMoleculesDeleted( DomainDecompBase* domainDecomposition){return _nNumMoleculesDeletedGlobalAlltime;}
	uint64_t getnNumMoleculesDeleted2( DomainDecompBase* domainDecomposition);

	void prepare_start(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer, double cutoffRadius);
	void init_positionMap(ParticleContainer* particleContainer);
	void preForce_action(ParticleContainer* particleContainer, double cutoffRadius);
	void postForce_action(ParticleContainer* particleContainer, DomainDecompBase* domainDecomposition);

	// connection to DensityControl
	void IncrementDeletedMoleculesLocal() {_nNumMoleculesDeletedLocal++;}
	void IncrementChangedMoleculesLocal() {_nNumMoleculesChangedLocal++;}
	void StoreDensity(const double& dVal) {_vecDensityValues.push_back(dVal);}
	void StoreValuesCV(const double& dDensity, const double& dVolume) {_dDensityTarget = dDensity; _dVolumeCV = dVolume;}

private:
	void ReadReservoir(DomainDecompBase* domainDecomp);
	void ReadReservoirFromFile(DomainDecompBase* domainDecomp);
	void ReadReservoirFromMemory(DomainDecompBase* domainDecomp);
	void DetermineMaxMoleculeIDs(DomainDecompBase* domainDecomp);
	void writeRestartfile();
	void calcDeltaY() { _dY = _dDeletedMolsPerTimestep * _dInvDensityArea; }
	void calcDeltaYbyDensity();
	// Check if molecule is a trapped one
	bool IsTrappedMolecule(const uint8_t& cid) {return cid != _vecChangeCompIDsUnfreeze.at(cid);}
	bool IsBehindTransitionPlane(const double& dPosY) {
		bool bRet = ( MD_LEFT_TO_RIGHT == _nMovingDirection && dPosY > _dTransitionPlanePosY ) ||
					( MD_RIGHT_TO_LEFT == _nMovingDirection && dPosY < _dTransitionPlanePosY );
		return bRet;
	}
	void InitSlabIndex() {
		switch(_nMovingDirection) {
		case MD_LEFT_TO_RIGHT:
			_nSlabindex = _reservoirSlabs-1; break;
		case MD_RIGHT_TO_LEFT:
			_nSlabindex = 0; break;
		}
	}
	void InitTransitionPlane(Domain* domain)
	{
		double dBoxY = domain->getGlobalLength(1);
		if(MD_LEFT_TO_RIGHT == _nMovingDirection)
			_dTransitionPlanePosY = 2*_dSlabWidth;
		else
			_dTransitionPlanePosY = dBoxY - 2*_dSlabWidth;
	}
	void NextReservoirSlab();
	void InsertReservoirSlab(ParticleContainer* particleContainer);

private:
	double _dDensityReservoir;
	double _dAreaXZ;
	double _dInvDensityArea;
	double _dY;
	double _dYInit;
	double _dYsum;
	double _velocityBarrier;
	double _dSlabWidthInit;
	double _dSlabWidth;
	double _dReservoirWidthY;
	uint64_t _nUpdateFreq;
	uint64_t _nWriteFreqRestart;
	uint64_t _nMaxMoleculeID;
	uint64_t _nMaxReservoirMoleculeID;
	uint64_t _nNumMoleculesDeletedLocal;
	uint64_t _nNumMoleculesDeletedGlobal;
	uint64_t _nNumMoleculesDeletedGlobalAlltime;
	uint64_t _nNumMoleculesChangedLocal;
	uint64_t _nNumMoleculesChangedGlobal;
	uint64_t _nNumMoleculesTooFast;
	uint64_t _nNumMoleculesTooFastGlobal;
	uint64_t _reservoirNumMolecules;
	uint64_t _reservoirSlabs;
	int32_t _nSlabindex;
	uint8_t _nReadReservoirMethod;
	uint8_t _nMovingDirection;
	uint8_t _nFeedRateMethod;
	std::string _reservoirFilename;
	std::map<uint64_t, std::array<double,10> > _storePosition;  //Map for frozen particle position storage <"id, position">
	std::vector< std::vector<Molecule> >_reservoir;
	bool _bIsRestart;  // simulation is a restart?
	std::list<uint64_t> _listDeletedMolecules;
	uint32_t _nNumValsSummation;
	uint64_t _numDeletedMolsSum;
	double _dDeletedMolsPerTimestep;
	double _dInvNumTimestepsSummation;
	bool _bMirrorActivated;
	double _dMirrorPosY;
	// identity change (by component ID)
	std::vector<uint32_t> _vecChangeCompIDsFreeze;
	std::vector<uint32_t> _vecChangeCompIDsUnfreeze;
	uint64_t _nDeleteNonVolatile;
	double _dMoleculeDiameter;
	double _dTransitionPlanePosY;
	// throttle parameters for each component
	std::vector<double> _vecThrottleFromPosY;
	std::vector<double> _vecThrottleToPosY;
	std::vector<double> _vecThrottleForceY;
	std::vector<double> _vecVeloctiyBarriers;
	// keep gas phase density
	std::vector<double> _vecDensityValues;
	double _dDensityTarget;
	double _dVolumeCV;
};

#endif /* METTDEAMON_H_ */
