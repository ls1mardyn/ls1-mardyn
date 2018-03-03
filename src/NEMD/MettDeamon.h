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
#include "utils/CommVar.h"

#include <map>
#include <array>
#include <fstream>
#include <vector>
#include <cstdint>
#include <limits>
#include <algorithm>

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
	FRM_DENSITY = 3,
	FRM_CONSTANT = 4
};

enum Zone2Method : uint8_t
{
	Z2M_UNKNOWN = 0,
	Z2M_RESET_ALL = 1,
	Z2M_RESET_YPOS_ONLY = 2,
};

enum MoleculeFormat : uint32_t {
	ICRVQD, IRV, ICRV
};

struct RestartInfoType
{
	uint32_t nBindindex;
	double dYsum;
};

struct FeedRateStruct
{
	uint32_t cid_target;

	struct FeedStruct
	{
		double init;
		double actual;
		double target;
		double sum;

	} feed;

	struct numMoleculesStruct
	{
		CommVar<uint64_t> inserted;
		CommVar<uint64_t> deleted;
		CommVar<uint64_t> changed_to;
		CommVar<uint64_t> changed_from;
	} numMolecules;
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
	uint32_t getFeedRateTargetComponentID() {return _feedrate.cid_target;}
	void IncrementInsertedMoleculesLocal() {_feedrate.numMolecules.inserted.local++;}
	void IncrementDeletedMoleculesLocal() {_feedrate.numMolecules.deleted.local++;}
	void IncrementChangedToMoleculesLocal() {_feedrate.numMolecules.changed_to.local++;}
	void IncrementChangedFromMoleculesLocal() {_feedrate.numMolecules.changed_from.local++;}
	void StoreDensity(const double& dVal) {_vecDensityValues.push_back(dVal);}
	void StoreValuesCV(const double& dDensity, const double& dVolume) {_dDensityTarget = dDensity; _dVolumeCV = dVolume;}

private:
	void findMaxMoleculeID(DomainDecompBase* domainDecomp);
	void writeRestartfile();
	void calcDeltaY();
	void calcDeltaYbyDensity();
	// Check if molecule is a trapped one
	bool IsTrappedMolecule(const uint8_t& cid) {return cid != _vecChangeCompIDsUnfreeze.at(cid);}
	bool IsBehindTransitionPlane(const double& dPosY) {
		bool bRet = ( MD_LEFT_TO_RIGHT == _nMovingDirection && dPosY > _dTransitionPlanePosY ) ||
					( MD_RIGHT_TO_LEFT == _nMovingDirection && dPosY < _dTransitionPlanePosY );
		return bRet;
	}
	bool IsInsideOuterReservoirSlab(const double& dPosY, const double& dBoxY);
	void releaseTrappedMolecule(Molecule* mol);
	void resetPositionAndOrientation(Molecule* mol, const double& dBoxY);

	void InitTransitionPlane(Domain* domain);
	void InsertReservoirSlab(ParticleContainer* particleContainer);
	void initRestart();

private:
	Reservoir* _reservoir;
	bool _bIsRestart;  // simulation is a restart?
	bool _bMirrorActivated;
	double _dAreaXZ;
	double _dInvDensityArea;
	double _velocityBarrier;
	double _dDeletedMolsPerTimestep;
	double _dInvNumTimestepsSummation;
	double _dMirrorPosY;
	double _dMoleculeDiameter;
	double _dTransitionPlanePosY;
	double _dDensityTarget;
	double _dVolumeCV;
	uint64_t _nUpdateFreq;
	uint64_t _nWriteFreqRestart;
	uint64_t _nMaxMoleculeID;
	uint64_t _nMaxReservoirMoleculeID;
	uint64_t _nNumMoleculesDeletedGlobalAlltime;
	uint64_t _nNumMoleculesTooFast;
	uint64_t _nNumMoleculesTooFastGlobal;
	uint8_t _nMovingDirection;
	uint8_t _nFeedRateMethod;
	uint8_t _nZone2Method;
	uint32_t _nNumValsSummation;
	int64_t _numDeletedMolsSum;
	uint64_t _nDeleteNonVolatile;
	std::map<uint64_t, std::array<double,10> > _storePosition;  //Map for frozen particle position storage <"id, position">
	std::list<uint64_t> _listDeletedMolecules;
	// identity change (by component ID)
	std::vector<uint32_t> _vecChangeCompIDsFreeze;
	std::vector<uint32_t> _vecChangeCompIDsUnfreeze;
	// keep gas phase density
	std::vector<double> _vecDensityValues;
	// throttle parameters for each component
	std::vector<double> _vecThrottleFromPosY;
	std::vector<double> _vecThrottleToPosY;
	std::vector<double> _vecThrottleForceY;
	std::vector<double> _vecVeloctiyBarriers;
	RestartInfoType _restartInfo;
	struct{
		double ymin;
		double ymax;
	} _manipfree;
	FeedRateStruct _feedrate;
};

struct FilepathStruct
{
	std::string data;
	std::string header;
};

struct BoxStruct
{
	double volume;
	std::array<double,3> length;
};

struct DensityStruct
{
	CommVar<uint64_t> numMolecules;
	double density;
};

class BinQueue;
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
	double getDensity(const uint32_t& cid) {return _density.at(cid).density;}
	void setDensity(const uint32_t& cid, const double& dVal) {_density.at(cid).density = dVal;}
	double getBoxLength(const uint32_t& nDim) {return _box.length.at(nDim);}
	void setBoxLength(const uint32_t& nDim, const double& dVal) {_box.length.at(nDim)=dVal;}
	double getVolume() {return _box.volume;}
	void setVolume(const double& dVal) {_box.volume = dVal;}
	double getBinWidth() {return _dBinWidth;}

	// queue methods
	uint32_t getActualBinIndex();
	uint64_t getNumMoleculesLocal();
	uint32_t getNumBins();
	std::vector<Molecule>& getParticlesActualBin();
	void nextBin(uint64_t& nMaxID);
	uint64_t getMaxMoleculeID();
	bool activateBin(uint32_t nBinIndex);

private:
	void calcPartialDensities(DomainDecompBase* domainDecomp);
	void changeComponentID(Molecule& mol, const uint32_t& cid);

private:
	MettDeamon* _parent;
	MoleculeDataReader* _moleculeDataReader;
	BinQueue* _binQueue;
	uint64_t _numMoleculesRead;
	uint64_t _nMaxMoleculeID;
	uint32_t _nMoleculeFormat;
	uint8_t _nReadMethod;
	double _dReadWidthY;
	double _dBinWidthInit;
	double _dBinWidth;
	std::vector<Molecule> _particleVector;
	std::vector<uint32_t> _vecChangeCompIDs;
	std::vector<DensityStruct> _density;
	FilepathStruct _filepath;
	BoxStruct _box;
};



class BinQueue
{
	class Bin
	{
	friend class BinQueue;
	public:
		Bin(std::vector<Molecule> vec, uint32_t nIndex) : _next(nullptr), _nIndex(nIndex)
		{
			for(auto p:vec)
				_particles.push_back(p);
		}
		Bin* _next;
		uint32_t _nIndex;
		std::vector<Molecule> _particles;
	};
private:
	Bin *_first, *_last, *_actual;
	uint32_t _numBins;
	uint32_t _nRoundCount;
	uint64_t _numParticles;
	uint64_t _maxID;


public:
	BinQueue() :
		_first(nullptr),
		_last(nullptr),
		_actual(nullptr),
		_numBins(0),
		_nRoundCount(0),
		_numParticles(0),
		_maxID(0)
	{
	}

	BinQueue(std::vector<Molecule> vec) :
		_first(nullptr),
		_last(nullptr),
		_actual(nullptr),
		_numBins(0),
		_nRoundCount(0),
		_numParticles(0),
		_maxID(0)
	{
		enque(vec);
	}

	~BinQueue() {
		_last->_next = nullptr;
		while (!isEmpty())
			deque();
	}

	bool isEmpty() {
		return _first == nullptr;
	}

	void enque(std::vector<Molecule> vec)
	{
		Bin* ptr = new Bin(vec, _numBins);
		if (isEmpty()) {
			_last = ptr;
			_first = ptr;
		} else {
			_last->_next = ptr;
			_last = ptr;
		}
		_numBins++;
		_numParticles += vec.size();
		// update max particle ID
		if(_numParticles > 0)
		{
			std::vector<Molecule>::iterator it = max_element(vec.begin(), vec.end(), molecule_id_compare);
			_maxID = ( (*it).id()>_maxID ? (*it).id() : _maxID);
		}
		else
			_maxID = 0;
		cout << "_maxID=" << _maxID << endl;
	}

	void deque()
	{
		if (!isEmpty()) {
			Bin* ptr = _first;
			_numParticles -= ptr->_particles.size();
			_first = _first->_next;
			delete ptr;
			_numBins--;
		}
		// update max particle ID
		this->determineMaxID();
	}

	std::vector<Molecule>& head() {
		return _first->_particles;
	}

	std::vector<Molecule>& getParticlesActualBin() {
		return _actual->_particles;
	}

	void showActualBin() {
		for(auto& p:_actual->_particles)
			std::cout << p << ", ";
		std::cout << endl;
	}

	void connectTailToHead()
	{
		_last->_next=_first;
		_actual = _first;
	}

	void next()
	{
		_actual = _actual->_next;
		if(_actual == _first)
			_nRoundCount++;
	}

	bool activateBin(uint32_t nBinIndex)
	{
		Bin* ptr = _first;
		while(ptr != nullptr)
		{
			cout << "ptr->_nIndex="<<ptr->_nIndex<<", nBinIndex="<<nBinIndex<<endl;
			if(ptr->_nIndex == nBinIndex)
				return true;
			ptr = ptr->_next;
			if(ptr == _first)
				break;
		}
		return false;
	}

	uint32_t getActualBinIndex() {return _actual->_nIndex;}
	uint32_t getNumBins() {return _numBins;}
	uint32_t getRoundCount() {return _nRoundCount;}
	uint64_t getNumParticles() {return _numParticles;}
	uint64_t getMaxID() {return _maxID;}

private:
	static bool molecule_id_compare(Molecule a, Molecule b)
	{
		return (a.id() < b.id());
	}
	void determineMaxID()
	{
		_maxID = 0;
		Bin* ptr = _first;
		while(ptr != nullptr)
		{
			std::vector<Molecule> vec = ptr->_particles;
			std::vector<Molecule>::iterator it = max_element(vec.begin(), vec.end(), molecule_id_compare);
			_maxID = ( (*it).id()>_maxID ? (*it).id() : _maxID);
			ptr = ptr->_next;
			if(ptr == _first)
				break;
		}
	}
};

#endif /* METTDEAMON_H_ */

