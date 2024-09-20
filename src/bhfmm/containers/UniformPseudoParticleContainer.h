/*
 * UniformPseudoParticleContainer.h
 *
 *  Created on: Feb 5, 2015
 *      Author: tchipevn
 */

#ifndef UNIFORMPSEUDOPARTICLECONTAINER_H_
#define UNIFORMPSEUDOPARTICLECONTAINER_H_

#include "PseudoParticleContainer.h"
#include "LeafNodesContainer.h"
#include "parallel/DomainDecompBase.h"
#include "bhfmm/utils/WignerMatrix.h"
#include "bhfmm/utils/RotationParameter.h"

#ifdef FMM_FFT
#include "bhfmm/fft/FFTAccelerationAPI.h"
#include "bhfmm/fft/FFTAccelerationAPI_extensions.h"
#include "bhfmm/fft/FFTSettings.h"
#include "bhfmm/fft/FFTFactory.h"
#include "bhfmm/fft/FFTOrderReduction.h"
#include "bhfmm/fft/TransferFunctionManagerAPI.h"
#include <bhfmm/FastMultipoleMethod.h>
#endif /* FMM_FFT */

#include <vector>
#include <map>
#include "bhfmm/HaloBufferNoOverlap.h"
#include "bhfmm/HaloBufferOverlap.h"

class Domain;
class DomainDecompBase;

namespace bhfmm {
class UniformPseudoParticleContainer: public PseudoParticleContainer {
public:
	UniformPseudoParticleContainer(double domainLength[3],
								   double bBoxMin[3],
								   double bBoxMax[3],
								   double LJCellLength[3],
								   unsigned LJSubdivisionFactor,
								   int orderOfExpansions,
								   ParticleContainer* ljContainer,
								   bool periodic = true
#ifdef QUICKSCHED
								   , qsched *scheduler = nullptr
#endif
	);
	~UniformPseudoParticleContainer();

	void clear();
	void build(ParticleContainer* pc);
	void upwardPass(P2MCellProcessor * cp);
	void horizontalPass(VectorizedChargeP2PCellProcessor * cp);
	void downwardPass(L2PCellProcessor *cp);

	// P2M
	void processMultipole(ParticleCellPointers& cell);

	// L2P
	void processFarField(ParticleCellPointers& cell);

	// M2M, M2L, L2L
	void processTree();

	//prints timer values to the standard output
	void printTimers();

    std::vector<std::vector<MpCell>> &getMpCellGlobalTop() ;

#ifdef FMM_FFT
    FFTAccelerationAPI *getFFTAcceleration() ;
#endif

    template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
    void M2LTowerPlateStep(int m1Loop, int mpCells, int curLevel);

// stuff used by Quicksched
    void M2MCompleteCell(int targetId, int level, int cellsPerDim);
    void P2MCompleteCell(int sourceId);
    void M2LCompleteCell(int targetId, int level, int cellsPerDimension);
    void M2LPair2Way(int cellA, int cellB, int level, int cellsPerDimension);
	void L2LCompleteCell(int sourceId, int level, int cellsPerDimension);
	void L2PCompleteCell(int targetId);
    enum taskModelTypesM2L {
        CompleteTarget,
        Pair2Way
    };

private:
	LeafNodesContainer* _leafContainer;
	int _wellSep;
	int _maxLevel;	//number of tree levels
	int _globalLevel;	//number of levels in global tree
	//StopLevel (only valid if avoidAllreduce=true):
	//level in global tree that decides where to change between global abd local allreduce
	//stopLevel = 1 -> only local reduces
	//stopLevel = globalLevel + 1 -> only global reduce
	//else local reduces from globalLevel till stopLevel and global reduce above
	int _stopLevel;
	//In the parallel version the octree is divided into two trees:
	//- A local subtree starting at _globalLevel + 1 which contains only the
	//ancestor nodes from the node that was assigned to the MPI process on
	//the _globalLevel
	//- A global tree which contains all the octree elements of the FMM tree
	//from level 0 to _globalLevel
	std::vector<std::vector<MpCell> > _mpCellGlobalTop;
	std::vector<std::vector<MpCell> > _mpCellLocal;
	double _cellLength[3];
	int _globalNumCellsPerDim;
	Domain* _domain;
	int _globalNumCells;	//total amount of cells
	std::vector<int> _occVector; // array for MPI allgather

	int _coeffVectorLength; //size of MPI buffer for multipole coefficients
	int _expansionSize; //size of one local or multipole expansion in doubles
	std::vector<double> _coeffVector; // array for MPI allgather
	std::vector<double> _coeffVector_me;
#ifdef ENABLE_MPI
	HaloBufferNoOverlap<double> * _multipoleRecBuffer, *_multipoleBuffer; //Buffer with use for non-overlapping communication
	HaloBufferOverlap<double> * _multipoleRecBufferOverlap, * _multipoleBufferOverlap, * _multipoleBufferOverlapGlobal, * _multipoleRecBufferOverlapGlobal; //Buffers for receiving and sending of global and local tree halos
	MPI_Request _allReduceRequest; //request that is used to check if Iallreduce of global reduce has finished
#endif
	bool _periodicBC;
	bool _avoidAllReduce; //if true then local reduces are performed according to stopping level; if false global allreduce is performed
	bool _importWholeGlobalRegion; //indicates if import for whole parent region should be imported (currently used when number of processor not a power of 2)
	bool _fuseGlobalCommunication; //indicates if fuse algorithm should be applied to reduce communication partners -> collect all values of parent region
	Vector3<int> _numProcessorsPerDim; //number of processors in every dimension
	Vector3<int> _numCellsOnGlobalLevel; //number of cells every processor owns on global level for each dimension (rectangular region)
	Vector3<int> _processorPositionGlobalLevel; //position of the processor on the global level
	Vector3<double> _bBoxMin; //minimum coordinates of the bounding box
	std::vector<int> _neighbours; //vector storing MPI ranks of neighbour ranks
	std::vector<std::vector<std::vector<int>>> _allRanks; //3 dimensional vector storing the whole 3 dimensional grid of mpi ranks
	int _globalLevelNumCells; //number of cells on the complete global level

#ifdef FMM_FFT
	TransferFunctionManagerAPI* _FFT_TM;
	FFTAccelerationAPI* _FFTAcceleration;
#endif  /* FMM_FFT */


	// M2M
	void CombineMpCell_Global(double *cellWid, int mpCells, int curLevel);

	// M2M
	void CombineMpCell_Local(double *cellWid, Vector3<int> localMpCells, int curLevel, Vector3<int> offset);

	// M2L
	void GatherWellSepLo_Global(double *cellWid, int mpCells, int curLevel);

	// M2L
	void GatherWellSepLo_Local(double *cellWid, Vector3<int> localMpCells, int curLevel, int doHalos);


#ifdef FMM_FFT
	// M2L
	void GatherWellSepLo_FFT_Global(double *cellWid, int mpCells, int curLevel);
	// M2L
	void GatherWellSepLo_FFT_Local(double *cellWid, Vector3<int> localMpCells, int curLevel, int doHalos);

	template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
	void GatherWellSepLo_FFT_Global_template(double *cellWid, int mpCells, int curLevel);

	template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
	void GatherWellSepLo_FFT_Local_template(double *cellWid, Vector3<int> localMpCells, int curLevel, int doHalos);
#endif /* FMM_FFT */

	// L2L
	void PropagateCellLo_Global(double *cellWid, int mpCells, int curLevel);

	// L2L
	void PropagateCellLo_Local(double *cellWid, Vector3<int> localMpCells, int curLevel, Vector3<int> offset);

	// for parallelization
	void AllReduceMultipoleMoments();
	void AllReduceLocalMoments(int mpCells, int _curLevel);
	/**
	 * Allreduce of all Multipole expansions from one level to the top is started
	 * Asynchronous collective communication
	 * @param mpCells: number of cells on curLevel
	 * @param curLevel: starting level
	 */
	void AllReduceMultipoleMomentsLevelToTop(int mpCells, int _curLevel);
	/**
	 * Allreduce of all Multipole expansions from one level to the top is completed by writing the received values
	 * @param mpCells: number of cells on curLevel
	 * @param curLevel: starting level
	 */
	void AllReduceMultipoleMomentsSetValues(int mpCells, int _curLevel);

	/**
	 * Reads multipole or local expansion values in a defined area into a buffer for each level in local tree (e.g. halo values)
	 * xLow,xHigh,yLow,yHigh,zLow and zHigh define the range of the rectangular subarea
	 * negative values address values at the end of the array (e.g. -1 means maxnumber - 1 in the respective dimension)
	 * 0 values refer to the beginning for low values or for the end of the array for the high values.
	 * @param doLocalExpansion: indicates if local or multipole expansions are read
	 * @param localMpCellsBottom: number of multipole cells per dimension at lowest level
	 * @param bottomLevel: lowest level
	 * @param buffer: buffer where values are should be stored
	 */
	void getHaloValues(Vector3<int> localMpCellsBottom,int bottomLevel, std::vector<double> buffer,
			int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh, bool doLocalExpansion);
	/**
	 * Sets multipole or local expansion values in a defined area into a buffer for each level in local tree (e.g. halo values)
	 * xLow,xHigh,yLow,yHigh,zLow and zHigh define the range of the rectangular subarea
	 * negative values address values at the end of the array (e.g. -1 means maxnumber - 1 in the respective dimension)
	 * 0 values refer to the beginning for low values or for the end of the array for the high values.
	 * @param doLocalExpansion: indicates if local or multipole expansions are written
	 * @param localMpCellsBottom: number of multipole cells per dimension at lowest level
	 * @param bottomLevel: lowest level
	 * @param buffer: buffer where values are stored in which should be written into the area
	 */
	void setHaloValues(Vector3<int> localMpCellsBottom,int bottomLevel, std::vector<double> bufferRec,
			int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh, bool doLocalExpansion);
	//for parallelization
	void communicateHalosNoOverlap();
	void communicateHalosOverlapStart(); //start communication with overlap for local tree
	void communicateHalosOverlapPostProcessingStart(); //backward communication start (asynchronous)
	void communicateHalosOverlapPostProcessingSetHalos(); //backward communication finish (set values)
	//has to be called after receive finished
	void communicateHalosOverlapSetHalos(); //finish communication with overlap (set values)
	void communicateHalos();
	/**
	 * Method used for sending global cell values that are calculated by MPI rank to other ranks that need to import them for M2L.
	 * Can also be used for receiving contributions of other MPI ranks to the own global values (receive = true in this case; used for NT method) (cells up the tree starting from cell on global level)
	 */
	void communicateOwnGlobalValue(int stopLevel = 1, bool receive = false);
	/**
	 * Method used for receiving halo values that are required for M2L calculation in global tree.
	 * Can also be used for sending results in halo cells to other MPI ranks to add them to their local values (send = true in this case; used in NT method) (cells up the tree starting from cell on global level)
	 */
	void communicateHaloGlobalValues(int stopLevel = 1, bool send = false);

	void communicateHalosX();
	void communicateHalosY();
	void communicateHalosZ();
	void communicateHalosAlongAxis(std::vector<double> lowerNeighbourBuffer, std::vector<double> higherNeighbourBuffer,
			std::vector<double> lowerNeighbourBufferRec, std::vector<double> higherNeighbourBufferRec,
			int lowerNeighbour, int higherNeighbour, int haloSize
			);
	bool _doNTLocal, _doNTGlobal; //indicate if NT method should be applied to the local tree part and/or the global tree part

	/**
	 * Performs a busy waiting scheme and return special integer values to indicate which MPI communications have already finished
	 * Needs to be called until all the communication in horizontal pass has finished (indicated by a return value of -1!
	 * For other definition of the other return values look into calling routine
	 */
	int busyWaiting();

	/**
	 * Initialize busy waiting variables
	 */
	void initBusyWaiting(){
		if((_avoidAllReduce && _stopLevel == 1) or _globalLevel == 0){
			_allReduceProcessed = 1;
		}
		else{
			_allReduceProcessed = 0;
		}
		_halosProcessed = 0;
		_sendLocalProcessed = 0;
		if(_globalLevel < 1 or !_avoidAllReduce or (_globalLevel == 1 and _fuseGlobalCommunication)){ // in this case only allreduce or nothing
			_sendGlobalProcessed = 1;
		}
		else{
			_sendGlobalProcessed = 0;
		}
		if(_doNTLocal){
			_backCommunicationLocalProcessed = 0;
		}
		else{
			_backCommunicationLocalProcessed = 1;
		}
		if(_doNTGlobal and _avoidAllReduce and not (_globalLevel == 1 and _fuseGlobalCommunication)){
			_backCommunicationGlobalProcessed = 0;
		}
		else{
			_backCommunicationGlobalProcessed = 1;
		}
		_backCommunicationLocalStarted = 0;
		_backCommunicationGlobalStarted = 0;
		if(_globalLevel < 1 or !_avoidAllReduce){
			_globalHalosProcessed = 1;
		}
		else{
			_globalHalosProcessed = 0;
		}
	}
	//filter methods that detect if certain combinations of cells are not allowed with current strategy (NT or not NT) for local and global tree parts
	bool filterM1Local(bool doHalos, int m1, int m1x, int m1y, int m1z, Vector3<int> localMpCells, int curLevel);
	bool filterM2Local(bool doHalos, int m1, int m1x, int m1y, int m1z, int m2, int m2x, int m2y, int m2z, Vector3<int> localMpCells, int curLevel, bool inHaloz, bool inHaloy, bool inHalox);
	bool filterM2Global(int curLevel, int *m2v, int *m1v, int m2x, int m2y, int m2z, int m2, int yOffset);

	//returns optimal value for the stopLevel of local to global Allreduce
	int optimizeAllReduce(/*ParticleContainer* ljContainer*/);
	//bool flags to indicate which communications have started and which have already been finished
	int _allReduceProcessed;
	int _globalHalosProcessed;
	int _halosProcessed; //local halos
	int _sendLocalProcessed;
	int _sendGlobalProcessed;
	int _backCommunicationLocalProcessed;
	int _backCommunicationLocalStarted;
	int _backCommunicationGlobalProcessed;
	int _backCommunicationGlobalStarted;

#ifdef ENABLE_MPI

	MPI_Comm _comm; //MPI communicator from domain decomposition
	MPI_Comm * _neighbourhoodComms; //MPI communicator that stores for each level the 8 MPI ranks that need to communicate for local reduces
	MPI_Comm * _allReduceComms; //MPI communicator that stores all MPI rank that need to communicate in global reduce for each possible stoplevel
#endif
	int _overlapComm; //indicates if overlap of communication is desired; Must be true currently!

#ifdef QUICKSCHED
    void generateResources(qsched *scheduler);
    /**
     * Needs M2L tasks
     */
    void generateP2MTasks(qsched* scheduler);
    /**
     * Needs M2L and P2M tasks
     */
    void generateM2MTasks(qsched* scheduler);
    /**
     * Needs L2P tasks
     */
    void generateP2PTasks(qsched *scheduler);
	/**
	 * Needs L2L tasks to be already created
	 */
    void generateM2LTasks(qsched *scheduler, taskModelTypesM2L taskModelM2L);
	/**
	 * Needs generateResources to be run first and L2P tasks created
	 */
    void generateL2LTasks(qsched *scheduler);
    /**
     * Need generateResources to be run first
     */
    void generateL2PTasks(qsched *scheduler);
#endif

};

} /* namespace bhfmm */

#endif /* UNIFORMPSEUDOPARTICLECONTAINER_H_ */
