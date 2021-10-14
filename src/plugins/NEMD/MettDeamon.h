/*
 * MettDeamon.h
 *
 *  Created on: 03.04.2017
 *      Author: thet
 */

#ifndef METTDEAMON_H_
#define METTDEAMON_H_

#include "../PluginBase.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "utils/CommVar.h"
//#include "NEMD/Request.h"

#include <map>
#include <array>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <memory>

// shuffle velocity list for each process before usage
#include <random>
#include <vector>
template < typename T > void shuffle( std::list<T>& lst ); // shuffle contents of a list

void create_rand_vec_ones(const uint64_t& nCount, const double& percent, std::vector<int>& v);
void update_velocity_vectors(std::unique_ptr<Random>& rnd, const uint64_t& numSamples, const double&T, const double&D, const double&v_neg, const double&e_neg,
		std::vector<double>& vxi, std::vector<double>& vyi, std::vector<double>& vzi);

#define FORMAT_SCI_MAX_DIGITS std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)

enum ReadReservoirMethods : uint8_t
{
	RRM_UNKNOWN = 0,
	RRM_READ_FROM_FILE = 1,
	RRM_READ_FROM_FILE_BINARY = 2,
	RRM_READ_FROM_MEMORY = 3,
	RRM_AMBIGUOUS = 4,
	RRM_EMPTY = 5,
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
	FRM_CONSTANT = 4,
	FRM_DIRECTED = 5
};

enum Zone2Method : uint8_t
{
	Z2M_UNKNOWN = 0,
	Z2M_RESET_ALL = 1,
	Z2M_RESET_YPOS_ONLY = 2,
};

enum ReleaseVelocityMethod : uint32_t
{
	RVM_UNKNOWN = 0,
	RVM_UNCHANGED = 1,
	RVM_FIX_VALUE = 2,
	RVM_ADD_FIX_VALUE = 3,
	RVM_NORM_DISTR = 4,
	RVM_NORM_DISTR_GENERATOR = 5
};

enum MoleculeFormat : uint32_t {
	ICRVQD, IRV, ICRV
};

struct RestartInfoType
{
	uint32_t nBindindex;
	double dYsum;
};

struct ParamsNormMB
{
	double temperature, drift;
	double a_neg, a_pos, v_neg, v_pos, e_neg, e_pos;
	std::string fpath;
};

struct FeedRateStruct
{
	uint32_t cid_target;
	uint32_t log_freq;

	struct FeedStruct
	{
		double init;
		double actual;
		double sum;

	} feed;

	struct ReleaseVelocityStruct
	{
		uint32_t method;
		double fix_value;
		ParamsNormMB normMB;
		std::vector<double> vx;
		std::vector<double> vy;
		std::vector<double> vz;
		uint64_t numSamples;
	} release_velo;

	struct numMoleculesStruct
	{
		CommVar<uint64_t> inserted;
		CommVar<uint64_t> deleted;
		CommVar<uint64_t> changed_to;
		CommVar<uint64_t> changed_from;
	} numMolecules;

	// rand vector, containing 100 values 0|1 for insertion part of trapped particles
	std::vector<int> vec_rand_ins;
};

class Domain;
class Ensemble;
class DomainDecompBase;
class ParticleContainer;
class XMLfileUnits;

class Reservoir;
class MettDeamon : public PluginBase
{
public:
	MettDeamon();
	~MettDeamon() override = default;

	/** @brief Read in XML configuration for MettDeamon and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
	<plugin name="MettDeamon">
		<control>
			<updatefreq>INT</updatefreq>                <!-- time span in which deleted particles are counted to determine feed rate -->
			<logfreqfeed>INT</logfreqfeed>             <!-- frequency with that the feed rate is logged -->
			<logfreqreleased>INT</logfreqreleased>     <!-- frequency with that velocity vectors of released particles are logged -->
			<writefreq>INT</writefreq>                <!-- write frequency of restart info (checkpoint) -->
			<numvals>INT</numvals>                       <!-- number of values collected to calculate an average feed rate -->
			<feed>
				<init>FLOAT</init>           <!-- initial feed rate  -->
				<direction>INT</direction>   <!-- 0: left --> right | 1: left <-- right  -->
				<method>INT</method>         <!-- feed rate method 4: fix rate | 5: get feedrate by MD Feedrate Director -->
				<targetID>INT</targetID>     <!-- component ID of particles feed rate determined from -->
				<target>FLOAT</target>       <!-- target value for feed rate, if method==4 -->
				<release_velo>
					<method>INT</method>           <!-- choose method for release velocity -->
					<fix_value>FLOAT</fix_value>   <!-- velocity with that particles are released from reservoir -->
				</release_velo>
			</feed>
			<z2method>1</z2method>                 <!-- choose zone2 method, 1:reset all i.e. also quaternion | 2:reset only y position od particles
			<manipfree> <ymin>50</ymin> <ymax>100</ymax> </manipfree>   <!-- range that is not affected with any manipulations -->
		</control>
		<reservoir update="1">   <!-- update="1": Update Reservoir's data structure before inserting new particles. This is mandatory when using kd-decomposition
			<file type="binary">
				<header>../../liq/run12/cp_binary-0.restart.header.xml</header>   <!-- checkpoint header file used for reservoir -->
				<data>../../liq/run12/cp_binary-0.restart.dat</data>              <!-- checkpoint data file used for reservoir -->
			</file>
			<binwidth>FLOAT</binwidth>         <!-- subdivision of reservoir phase into bins, overall reservoir width = 2*binwidth -->
			<ins_percent>FLOAT</ins_percent>   <!-- only the fraction ins_percent=FLOAT of reservoir particles will be inserted -->
		</reservoir>
		<restart>
			<binindex>INT</binindex>   <!-- index of last inserted (actual) reservoir bin  -->
			<deltaY>FLOAT</deltaY>     <!-- last inserted (actual) reservoir bin was moved forward already by deltaY=FLOAT -->
		</restart>

		<changes>
			<change> <from>INT</from> <to>INT</to> </change>   <!-- change component ID of reservoir particles from (moving) INT to (freeze) INT -->
		</changes>
	</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;
	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
	void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
	void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
	void siteWiseForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override {};
	void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override {};
	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {};
	std::string getPluginName() override {return std::string("MettDeamon");}
	static PluginBase* createInstance() {return new MettDeamon();}

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

	// connection to other general plugins
	void setActualFeedrate(const double& feed_actual) {
		if (FRM_DIRECTED == _nFeedRateMethod) {
			_feedrate.feed.actual = feed_actual;
			global_log->info() << "[MettDeamon]: Set new feed rate by MDFRD to vf= " << _feedrate.feed.actual << std::endl;
		} else {
			global_log->warning() << "[MettDeamon]: Feed rate not set because feed method ( " << _nFeedRateMethod << " ) is not " << FRM_DIRECTED << std::endl;
		}
	}
	void setInitFeedrate(const double& feed_init) {
		if (FRM_DIRECTED == _nFeedRateMethod) {
			_feedrate.feed.init = feed_init;
			global_log->info() << "[MettDeamon]: Set init feed rate by MDFRD to vf= " << _feedrate.feed.init << std::endl;
		} else {
			global_log->warning() << "[MettDeamon]: Feed rate not set because feed method ( " << _nFeedRateMethod << " ) is not " << FRM_DIRECTED << std::endl;
		}
	}
	double getInvDensityArea() {return _dInvDensityArea;}

private:
	void findMaxMoleculeID(DomainDecompBase* domainDecomp);
	void writeRestartfile();
	void logFeedrate();
	void logReleased();
	void logReleasedVelocities();
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
	void releaseTrappedMolecule(Molecule* mol, bool& bDeleteParticle);
	void resetPositionAndOrientation(Molecule* mol, const double& dBoxY);
	void resetVelocity(Molecule* mol);

	void InitTransitionPlane(Domain* domain);
	void getAvailableParticleIDs(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			CommVar<std::vector<uint64_t> >& particleIDs_available, const CommVar<uint64_t>& numParticleIDs);
	void updateReservoir(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer);
	void InsertReservoirSlab(ParticleContainer* particleContainer);
	void initRestart();
	
	// stat. evap
	void readNormDistr();

	// rand vector for trapped particle insertion
	void updateRandVecTrappedIns();

private:
	std::unique_ptr<Reservoir> _reservoir;
	bool _bIsRestart;  // simulation is a restart?
	bool _bInitFeedrateLog;
	bool _bInitRestartLog;
	double _dAreaXZ;
	double _dInvDensityArea;
	double _dDeletedMolsPerTimestep;
	double _dInvNumTimestepsSummation;
	double _dTransitionPlanePosY;
	double _dDensityTarget;
	double _dVolumeCV;
	uint64_t _nUpdateFreq;
	uint64_t _nWriteFreqRestart;
	CommVar<uint64_t> _nMaxMoleculeID;
	uint64_t _nMaxReservoirMoleculeID;
	uint64_t _nNumMoleculesDeletedGlobalAlltime;
	CommVar<uint64_t> _nNumMoleculesTooFast;
	uint8_t _nMovingDirection;
	FeedRateMethod _nFeedRateMethod;
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

	RestartInfoType _restartInfo;
	struct{
		double ymin, ymax;
		uint32_t cid_ub;
	} _manipfree;
	FeedRateStruct _feedrate;
	
	// stat. evap.
	struct{
		bool enabled;
		uint32_t type;  // 1|2 : left|right or liq|vap
		double pos_y;
		struct{
			struct{
				uint32_t pos;
				uint32_t neg;
			} left;
			struct{
				uint32_t pos;
				uint32_t neg;
			} right;
		} cid;
	} _vap_trans_plane;
	
	struct NormMB{
		struct NormFnames{
			std::string vxz;
			std::string vy;
		} fname;
		std::list<double> vxz;
		std::list<double> vy;
	} _norm;
	
	struct{
		CommVar<uint64_t> count;
		CommVar<uint64_t> deleted;
		uint64_t log_freq;
		uint64_t log_freq_vel;
		bool init_file;
		bool init_file_vel;
		std::vector<std::array<double,3> > log_v;
	} _released;

	std::unique_ptr<Random> _rnd;
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
	~Reservoir();

	void readXML(XMLfileUnits& xmlconfig);

	// read particle data
	void readParticleData(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer);
	void updateParticleData(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer);
private:
	void readFromMemory(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer);
	void readFromFile(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer);
	void readFromFileBinary(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer);
	void readFromFileBinaryHeader();
	void sortParticlesToBins(DomainDecompBase* domainDecomp, ParticleContainer* particleContainer);

public:
	// Getters, Setters
	double getDensity(const uint32_t& cid) {return _density.at(cid).density;}
	void setDensity(const uint32_t& cid, const double& dVal) {_density.at(cid).density = dVal;}
	double getBoxLength(const uint32_t& nDim) {return _box.length.at(nDim);}
	void setBoxLength(const uint32_t& nDim, const double& dVal) {_box.length.at(nDim)=dVal;}
	double getVolume() {return _box.volume;}
	void setVolume(const double& dVal) {_box.volume = dVal;}
	double getBinWidth() {return _dBinWidth;}
	double GetInsPercent() {return _dInsPercent;}
	void setInsPercent(const double& dVal) {_dInsPercent = dVal;}
	uint8_t getReadMethod() {return _nReadMethod;}

	// queue methods
	uint32_t getActualBinIndex();
	uint64_t getNumMoleculesLocal();
	uint32_t getNumBins();
	std::vector<Molecule>& getParticlesActualBin();
	bool nextBin(uint64_t& nMaxID);
	uint64_t getMaxMoleculeID();
	bool activateBin(uint32_t nBinIndex);
	void clearBinQueue();
	void printBinQueueInfo();

private:
	void calcPartialDensities(DomainDecompBase* domainDecomp);
	void changeComponentID(Molecule& mol, const uint32_t& cid);
	bool isRelevant(DomainDecompBase* domainDecomp, Domain* domain, Molecule& mol);

private:
	MettDeamon* _parent;
	std::unique_ptr<MoleculeDataReader> _moleculeDataReader;
	std::unique_ptr<BinQueue> _binQueue;
	uint64_t _numMoleculesRead;
	uint64_t _nMaxMoleculeID;
	uint32_t _nMoleculeFormat;
	uint8_t _nReadMethod;
	double _dReadWidthY;
	double _dBinWidthInit;
	double _dBinWidth;
	double _dInsPercent;  // only insert percentage of reservoir density
	std::vector<Molecule> _particleVector;
	std::vector<uint32_t> _vecChangeCompIDs;
	std::vector<DensityStruct> _density;
	FilepathStruct _filepath;
	BoxStruct _box;
	bool _bUpdateBinQueue;  // BinQueue have to be updated if bounding boxes changes, e.g. in case of using kd-decomposition
};



class BinQueue
{
	class Bin
	{
	friend class BinQueue;
	public:
		Bin(const std::vector<Molecule>& vec, uint32_t nIndex) : _next(nullptr), _nIndex(nIndex)
		{
			for(const auto& p:vec)
				_particles.push_back(p);
		}
		uint32_t getIndex() {return _nIndex;}
		uint64_t getNumParticles() {return _particles.size();}
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
		enque(std::move(vec));
	}

	~BinQueue() {
		if(_last){
			// break connection of last element to first element (established via connectTailToHead())
			_last->_next = nullptr;
		}
		auto ptr = _first;
		while (ptr) {
			// as long as ptr isn't nullptr we continue to delete the next element.
			auto ptr_next = ptr->_next;
			delete ptr;
			ptr = ptr_next;
		}
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
			_actual = ptr;
		} else {
			_last->_next = ptr;
			_last = ptr;
		}
		_last->_next = _first;  // connect tail to head
		_numBins++;
		_numParticles += vec.size();
		// update max particle ID
		if(_numParticles > 0)
		{
			std::vector<Molecule>::iterator it;
			if(vec.empty())
				it = vec.end();
			else if(vec.size() == 1)
				it = vec.begin();
			else
				it = max_element(vec.begin(), vec.end(), molecule_id_compare);

			if(it != vec.end() )
				_maxID = ( (*it).getID()>_maxID ? (*it).getID() : _maxID);
		}
		else
			_maxID = 0;
	}

	void deque()
	{
		if (isEmpty()) {
			return;
		}
		else if (_first == _last) {
			_numParticles -= _first->getNumParticles();
			delete _first;
			_first = _last = nullptr;
			_numBins--;
		}
		else {
			Bin* ptr = _first;
			while(ptr->_next != _last) {
				ptr = ptr->_next;
			}
			_numParticles -= ptr->_next->getNumParticles();
			delete ptr->_next;
			ptr->_next = nullptr;
			_last = ptr;
			_numBins--;
			_last->_next = _first;  // connect tail to head
		}
	}
	
	void clear()
	{
		while (!isEmpty()) {
			deque();
		}
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

	bool next()
	{
		_actual = _actual->_next;
		if(_actual == _first)
			_nRoundCount++;
		bool bSuccess = dynamic_cast<Bin*>(_actual);
		return bSuccess;
	}

	bool activateBin(uint32_t nBinIndex)
	{
		Bin* ptr = _first;
		while(ptr != nullptr)
		{
			if(ptr->_nIndex == nBinIndex) {
				_actual = ptr;
				return true;
			}
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
		return (a.getID() < b.getID());
	}
	void determineMaxID()
	{
		_maxID = 0;
		Bin* ptr = _first;
		while(ptr != nullptr)
		{
			std::vector<Molecule> vec = ptr->_particles;
			auto it = max_element(vec.begin(), vec.end(), molecule_id_compare);
			_maxID = ( (*it).getID()>_maxID ? (*it).getID() : _maxID);
			ptr = ptr->_next;
			if(ptr == _first)
				break;
		}
	}
};

#endif /* METTDEAMON_H_ */

