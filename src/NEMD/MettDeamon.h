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
	RRM_READ_FROM_FILE_BINARY = 2,
	RRM_READ_FROM_MEMORY = 3,
	RRM_AMBIGUOUS = 4,
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

enum MoleculeFormat : uint32_t {
	ICRVQD, IRV, ICRV
};

class Domain;
class Ensemble;
class DomainDecompBase;
class ParticleContainer;
class XMLfileUnits;

class Reservoir;
class MettDeamon
{
public:
	MettDeamon();
	~MettDeamon();

	void readXML(XMLfileUnits& xmlconfig);

	uint64_t getnNumMoleculesDeleted( DomainDecompBase* domainDecomposition){return _nNumMoleculesDeletedGlobalAlltime;}
	uint64_t getnNumMoleculesDeleted2( DomainDecompBase* domainDecomposition);
	uint8_t getMovingDirection() {return _nMovingDirection;}
	double  getTransitionPlanePosY() {return _dTransitionPlanePosY;}

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
	void findMaxMoleculeID(DomainDecompBase* domainDecomp);
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

	void InitTransitionPlane(Domain* domain);
	void InsertReservoirSlab(ParticleContainer* particleContainer);

private:
	double _dAreaXZ;
	double _dInvDensityArea;
	double _dY;
	double _dYInit;
	double _dYsum;
	double _velocityBarrier;
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
	uint8_t _nMovingDirection;
	uint8_t _nFeedRateMethod;
	std::map<uint64_t, std::array<double,10> > _storePosition;  //Map for frozen particle position storage <"id, position">
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
	Reservoir* _reservoir;
};

class MoleculeDataReader;
class Reservoir
{
public:
	Reservoir(MettDeamon* parent);
	~Reservoir(){}

	void readXML(XMLfileUnits& xmlconfig);

	// read particle data
	void readParticleData(DomainDecompBase* domainDecomp);
private:
	void readFromMemory(DomainDecompBase* domainDecomp);
	void readFromFile(DomainDecompBase* domainDecomp);
	void readFromFileBinary(DomainDecompBase* domainDecomp);
	void readFromFileBinaryHeader();
	void sortParticlesToBins();

public:
	// Getters, Setters
	uint64_t getNumMoleculesLocal() {return _numMoleculesLocal;}
	uint64_t getNumMoleculesGlobal() {return _numMoleculesGlobal;}
//	void setNumMolecules(uint64_t nVal) {_numMolecules = nVal;}
	uint64_t getNumBins() {return _numBins;}
	void setNumBins(uint32_t nVal) {_numBins = nVal; _binVector.resize(nVal);}
	std::string getFilename() {return _strFilename;}
	std::string getFilenameHeader() {return _strFilenameHeader;}
//	void setNumMolecules(uint64_t nVal) {_numMolecules = nVal;}
	double getDensity() {return _dDensity;}
	void setDensity(double dVal) {_dDensity = dVal;}
	double getBoxLength(uint32_t nDim) {return _arrBoxLength.at(nDim);}
	void setBoxLength(uint32_t nDim, double dVal) {_arrBoxLength.at(nDim)=dVal;}
	double getVolume() {return _dVolume;}
	void setVolume(double dVal) {_dVolume = dVal;}
	std::vector<Molecule>& getBinMoleculeVector(uint32_t nBinIndex) {return _binVector.at(nBinIndex);}
	std::vector<Molecule>& getBinMoleculeVectorActual() {return _binVector.at(_nBinIndex);}
	int32_t getBinIndex() {return _nBinIndex;}
	void setBinIndex(int32_t nVal) {_nBinIndex = nVal;}
	double getBinWidth() {return _dBinWidth;}

	// more methods
	void initBinIndex(uint8_t nMovingDirection);
	void nextBin(uint8_t nMovingDirection, uint64_t& nMaxID);
	uint64_t findMaxMoleculeID();

private:
	uint64_t calcNumMoleculesLocal();
	uint64_t calcNumMoleculesGlobal(DomainDecompBase* domainDecomp);

private:
	MettDeamon* _parent;
	uint64_t _numMoleculesRead;
	uint64_t _numMoleculesLocal;
	uint64_t _numMoleculesGlobal;
	uint64_t _numBins;
	uint64_t _nMaxMoleculeID;
	uint32_t _nMoleculeFormat;
	uint8_t _nReadMethod;
	int32_t _nBinIndex;
	double _dReadWidthY;
	double _dBinWidthInit;
	double _dBinWidth;
	double _dDensity;
	double _dVolume;
	std::string _strFilename;
	std::string _strFilenameHeader;
	MoleculeDataReader* _moleculeDataReader;
	std::array<double,3> _arrBoxLength;
	std::vector<Molecule> _particleVector;
	std::vector< std::vector<Molecule> > _binVector;
};

#endif /* METTDEAMON_H_ */
