#ifndef KDDECOMPOSITION2_H_
#define KDDECOMPOSITION2_H_

#include <mpi.h>
#include <algorithm>
#include <list>
#include <vector>

#define KDDIM 3

#include "DomainDecompMPIBase.h"
#include "parallel/CollectiveCommunication.h"

class ParticleData;
class KDNode;


/** @brief KD tree based domain decomposition for better load balancing.
 *
 * The basic idea is to build up all possible subdivisions and do a A*-like
 * search of the best subdivision.
 *
 * The function to minimize is the load imbalance:
 * sum over all children i: \sum (load_i - _optimalLoad)^2
 *
 * During the downward pass an estimate is the expected load imbalance, which is
 * averaged over all children, in the upward pass, the exact imbalance is calculated.
 *
 * \note It is important for the A*-search that the estimate is an
 *       underestimation (<=) of the load imbalance.
 *
 * \note Some computation of the deviation / expected deviation is done in KDNode.
 */
class KDDecomposition: public DomainDecompMPIBase {

	friend class KDDecompositionTest;

 public:
	/** @brief create an initial decomposition tree
	 *
	 * The constructor determines the number of cells and creates an initial decomposition
	 * of the domain (not yet balanced), which is stored in _decompTree and _ownArea.
	 * @param cutoffRadius largest cutoff radius of a molecule (determines a basic
	 *                     cell for loadbalancing)
	 * @param domain
	 * @param updateFrequency every n-th timestep, load will be balanced.
	 * @param fullSearchThreshold If a KDNode has a processor count less or equal this number,
	 *                            all possible decompositions will be investigated, so it
	 *                            influences the quality of the load balancing. I recommend to
	 *                            set it to 2 - 4.
	 */
	KDDecomposition(double cutoffRadius, Domain* domain, int updateFrequency = 100, int fullSearchThreshold = 2, bool hetero=false, bool cutsmaller=false, bool forceRatio=false);

	KDDecomposition();

	// documentation see father class (DomainDecompBase.h)
	~KDDecomposition();


	/** @brief Read in XML configuration for KDDecomposition and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <parallelisation type="KDDecomposition">
	     <updateFrequency>INTEGER</updateFrequency>
	     <fullSearchThreshold>INTEGER</fullSearchThreshold>
	   </parallelisation>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig) override;

	//###############################################
	//### The following methods are those of the  ###
	//### base class which have to be implemented ###
	//###############################################

	// documentation in base class
	virtual int getNonBlockingStageCount() override;

	// documentation in base class
	virtual void prepareNonBlockingStage(bool forceRebalancing,
				ParticleContainer* moleculeContainer, Domain* domain,
				unsigned int stageNumber) override;

	// documentation in base class
	virtual void finishNonBlockingStage(bool forceRebalancing,
			ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber) override;

	// documentation in base class
	bool queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain) override;

	void balanceAndExchange(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain) override;

	//! @todo comment and thing
	double getBoundingBoxMin(int dimension, Domain* domain) override;
	//! @todo comment and thing
	double getBoundingBoxMax(int dimension, Domain* domain) override;

	//! @brief writes information about the current decomposition into the given file
	//!
	//! This decomposition first writes a very small header with the following lines:
	//! - "size", followed by three double values for the size of the domain
	//! - "decompData Regions": just this text to state that now the header is finished and data follows
	//! - one line per proc, giving the bounding box (minx, miny, minz, maxx, maxy, maxz)
	//!
	//! An example file for 5 procs could look like this:
	//!
	//! |size 62.0 62.0 62.0 \n
	//!  decompData Regions \n
	//!  0.0 0.0 0.0 20.0 62.0 62.0 \n
	//!  20.0 0.0 0.0 62.0 40.0 25.0 \n
	//!  20.0 40.0 0.0 62.0 62.0 25.0 \n
	//!  20.0 0.0 25.0 62.0 30.0 62.0 \n
	//!  20.0 30.0 25.0 62.0 62.0 62.0
	//! @param filename name of the file into which the data will be written
	//! @param domain e.g. needed to get the bounding boxes
	void printDecomp(std::string filename, Domain* domain) override;

	int getUpdateFrequency() { return _frequency; }
	void setUpdateFrequency(int frequency) { _frequency = frequency; }
	virtual std::vector<int> getNeighbourRanks() override {
		//global_log->error() << "not implemented \n";
		exit(-1);
		return std::vector<int> (0);
	}
	virtual std::vector<int> getNeighbourRanksFullShell() override{
		//global_log->error() << "not implemented \n";
		exit(-1);
		return std::vector<int> (0);
	}

	virtual std::vector<CommunicationPartner> getNeighboursFromHaloRegion(Domain* domain, const HaloRegion& haloRegion, double cutoff) override;

 private:
	void constructNewTree(KDNode *& newRoot, KDNode *& newOwnLeaf, ParticleContainer* moleculeContainer);
	/**
	 *
	 * @param newRoot
	 * @param newOwnLeaf
	 * @param moleculeContainer
	 * @return true if OK, false if deadlock
	 */
	bool migrateParticles(const KDNode& newRoot, const KDNode& newOwnLeaf, ParticleContainer* moleculeContainer) const;
	void initCommunicationPartners(double cutoffRadius, Domain * domain);


	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//$ sonstige Methoden
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	void calculateCostsPar(KDNode* area, std::vector<std::vector<double> >& costsLeft, std::vector<std::vector<double> >& costsRight, MPI_Comm commGroup);


	//! @brief calculates the index of a certain cell in the global cell array
	//!
	//! During the construction of a new decomposition, a area is recursively
	//! divided into two smaller areas. Each of the corresponding dividing planes
	//! can be orthogonal to each of the three spacial axis. For each of these three
	//! divisions, all cells in the "left" part are traversed by using three loops over
	//! the three dimensions. E.g. if the dividing plane is orthogonal to the x-axis,
	//! the loop goes over the left part of all x-cells, over all y-cells and over all z-cells.
	//! To avoid code duplication for the three possibilities of the dividing plane, a
	//! special method is used. The tree loops don't run over x, y and z, but over
	//! divDim, dim1 and dim2, where those three values are chosen differently
	//! for each of the dividing planes (x: 0,1,2; y: 1,0,2; z: 2,0,1). Those three values plus
	//! the corresponding loop-indeces are used by this method to calculate the global index.
	unsigned int getGlobalIndex(int divDim, int dim1, int dim2, int index_dim, int index_dim1, int index_dim2, KDNode* localArea);

	//! @brief provides a "real" modulo function
	//!
	//! The C build-in modulo operation uses a symmetric modulo,
	//! This method provides a mathematical modulo (important for negative numbers)
	//! @param number Number to calculate the modulo for
	//! @param modulo obviously the modulo
	int ownMod(int number, int modulo) const;

	void getCellBorderFromIntCoords(double * lC, double * hC, int lo[3], int hi[3]) const;

	void getCellIntCoordsFromRegionPeriodic(int * lo, int * hi, const double lC[3], const double hC[3], const Domain* domain) const;

	//! @brief exchange decomposition data and build up the resulting tree
	//!
	//! After the new decomposition has been determined (by recDecompPar), each
	//! process knows its own Area (and that part of the decomp tree which lies on the
	//! the path between the root node and the node with the own Area), but each process
	//! has to know the complete decomposition. This method exchanges the decomposition
	//! data between the processes and builds up a complete decomposition Tree on
	//! all processes.
	void completeTreeInfo(KDNode*& root, KDNode*& ownArea);


	//! @brief Find all processes which own a part of the given area
	//!
	//! This method is supposed to run recursively. It should usually be called
	//! with decompTree and testNode pointing both to the root of a decomposition tree.
	//! In the recursive calls, testNode will be changed to walk through the whole tree.
	//! The given area can be partly outside of the area covered by the decomposition tree
	//! (usually because of the halo region), in this case the periodic boundary conditions
	//! are applied to get the corresponding area on the other side of the region.
	//! @param low low corner of the area for which the owners shall be found
	//! @param high high corner of the area for which the owners shall be found
	//! @param decompTree root KDNode of a domain decomposition tree. The area covered by
	//!                   this KDNode must not overlap the domain
	//! @param testNode KDNode currenty tested for intersection (method is used recursively)
	//! @param procIDs here the ids of those procs which own a part will be stored
	//! @param neighbHaloAreas here the corresponding (see procIDs) areas will be stored. For
	//!                        each procID six double values are stored, first the three
	//!                        coordinates for the low corner of the area, then the high corner
	//! @todo make it work for overlapping decomposition trees
	void getOwningProcs(int low[KDDIM], int high[KDDIM], KDNode* decompTree, KDNode* testNode, std::vector<int>* procIDs, std::vector<int>* neighbHaloAreas) const;


	//! @brief
	//! @todo _numParticles should perhaps not be a member variable (think about that)
	void getNumParticles(ParticleContainer* moleculeContainer);

	bool decompose(KDNode* fatherNode, KDNode*& ownArea, MPI_Comm commGroup);

	bool decompose(KDNode* fatherNode, KDNode*& ownArea, MPI_Comm commGroup, const double globalMinimalDeviation);

	/**
	 * Get the load for the area represented by this node.
	 *
	 * @note: if this implementation is too slow for large domains, we could change
	 *        it so that all division costs for
	 */
	bool calculateAllSubdivisions(KDNode* node, std::list<KDNode*>& subdivededNodes, MPI_Comm commGroup);
	
	
	void updateMeanProcessorSpeeds(std::vector<double>& processorSpeeds,
			std::vector<double>& accumulatedProcessorSpeeds, ParticleContainer* moleculeContainer);

	//######################################
	//###    private member variables    ###
	//######################################


	//! Length of one cell. The length of the cells is the cutoff-radius, or a
	//! little bit larger as there has to be a natural number of cells.
	double _cellSize[KDDIM];

	//! Number of cells in the global domain in each coordinate direction.
	//!  _globalCellsPerDim does not include the halo region
	int _globalCellsPerDim[KDDIM];

	//! global number of cells (without halo)
	int _globalNumCells;

	// Pointer to the root element of the decomposition tree
	KDNode* _decompTree;

	// each process owns an area in the decomposition
	KDNode* _ownArea;

	//! Number of particles for each cell (including halo?)
	unsigned int* _numParticlesPerCell;

	/* TODO: This may not be equal to the number simulation steps if balanceAndExchange
	 * is not called exactly once in every simulation step! */
	//! number of simulation steps. Can be used to trigger load-balancing every _frequency steps
	size_t _steps;

	//! determines how often rebalancing is done
	int _frequency;

	double _cutoffRadius;

	/*
	 * Threshold for full tree search. If a node has more than _fullSearchThreshold processors,
	 * it is for each dimension divided in the middle only. Otherwise, all possible subdivisions
	 * are created. // TODO hetero: change description
	 */
	int _fullSearchThreshold;


	std::vector<double> _processorSpeeds;
	std::vector<double> _accumulatedProcessorSpeeds;//length nprocs+1, first element is 0.
	double _totalMeanProcessorSpeed;
	double _totalProcessorSpeed;
	int _processorSpeedUpdateCount;
	bool _heterogeneousSystems;
	bool _splitBiggest;  // indicates, whether a subdomain is to be split along its biggest size
	bool _forceRatio;  // if you want to enable forcing the above ratio, enable this.
};


#endif /* KDDECOMPOSITION2_H_ */
