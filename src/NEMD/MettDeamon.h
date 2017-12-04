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

enum MoleculeFormat : uint32_t {
	ICRVQD, IRV, ICRV
};

struct RestartInfoType
{
	uint32_t nBindindex;
	double dYsum;
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
	void initRestart();

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
	RestartInfoType _restartInfo;
	double _dFeedRate;
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
	uint64_t getNumMoleculesGlobal() {return _numMoleculesGlobal;}
//	void setNumMolecules(uint64_t nVal) {_numMolecules = nVal;}
	std::string getFilename() {return _strFilename;}
	std::string getFilenameHeader() {return _strFilenameHeader;}
//	void setNumMolecules(uint64_t nVal) {_numMolecules = nVal;}
	double getDensity() {return _dDensity;}
	void setDensity(double dVal) {_dDensity = dVal;}
	double getBoxLength(uint32_t nDim) {return _arrBoxLength.at(nDim);}
	void setBoxLength(uint32_t nDim, double dVal) {_arrBoxLength.at(nDim)=dVal;}
	double getVolume() {return _dVolume;}
	void setVolume(double dVal) {_dVolume = dVal;}
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
	uint64_t calcNumMoleculesGlobal(DomainDecompBase* domainDecomp);

private:
	MettDeamon* _parent;
	uint64_t _numMoleculesRead;
	uint64_t _numMoleculesGlobal;
	uint64_t _nMaxMoleculeID;
	uint32_t _nMoleculeFormat;
	uint8_t _nReadMethod;
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
	BinQueue* _binQueue;
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
		uint32_t _nIndex;
		Bin* _next;
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
		std::vector<Molecule>::iterator it = max_element(vec.begin(), vec.end(), molecule_id_compare);
		_maxID = ( (*it).id()>_maxID ? (*it).id() : _maxID);
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

