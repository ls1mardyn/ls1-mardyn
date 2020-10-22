#ifndef KDNODE_H_
#define KDNODE_H_

// include because of macro KDDIM
#include "parallel/KDDecomposition.h" /* TODO seriously..? */

#include "parallel/MPIKDNode.h"

#include <vector>

class VTKGridWriterImplementation;

//! @brief represents a node in the decomposition tree when using KDDecomposition
//! @author Martin Buchholz, Wolfgang Eckhardt
//! 
//! The KDDecomposition decomposes the domain by recursively splitting the domain
//! into smaller parts. This class is used to represent this decomposition.
//! The root node of the decomposition covers the whole domain, in the first
//! splitting step, this domain is divided into two parts, where each part
//! then has to be divided into several smaller parts. How many parts/regions there
//! depends on the number of processes, each process will get one region.
//! So the KNNode also has to store how many process "share" the current region
//! The leaf nodes of the tree represent the region of the single processes.
//! The regions (the size of the regions) is not stored in some floating-point length unit
//! but in cells. So it is assumed that the domain is discretised with cells and
//! the decomposition is based on distributing those cells (blocks of cells) to the
//! processes
class KDNode {

public:

	/**
	 * MPIKDNodes represent the data of a KDNode which has to be sent via MPI.
	 * MPIKDNodes can be used directly in MPI Send/Receive operations, with
	 * mpi_data_type as type.
	 */
	typedef MPIKDNodePacked MPIKDNode;

	KDNode() :_numProcs(0), _nodeID(0), _owningProc(0),
			  _child1(nullptr), _child2(nullptr), _load(0.0), _optimalLoadPerProcess(0.0),
			  _expectedDeviation(0.0), _deviation(0.0), _level(0)
	{
	}


	/**
	 * Copy constructor copies everything except for the children (are set to nullptr!)
	 */
	KDNode(const KDNode& other) : _numProcs(other._numProcs), _nodeID(other._nodeID),
			_owningProc(other._owningProc), _child1(nullptr), _child2(nullptr),
			_load(other._load), _optimalLoadPerProcess(other._optimalLoadPerProcess),
			_expectedDeviation(other._expectedDeviation), _deviation(other._deviation),
			_level(other._level)
	{
		for (int dim = 0; dim < KDDIM; dim++) {
			_lowCorner[dim] = other._lowCorner[dim];
			_highCorner[dim] = other._highCorner[dim];
			_coversWholeDomain[dim] = other._coversWholeDomain[dim];
		}
	}

	KDNode(int numP, const int low[KDDIM], const int high[KDDIM], int id, int owner, bool coversAll[KDDIM], int level)
	: _numProcs(numP), _nodeID(id), _owningProc(owner),
	  _child1(nullptr), _child2(nullptr), _load(0.0), _optimalLoadPerProcess(0.0),
	  _expectedDeviation(0.0), _deviation(0.0), _level(level)
	{
		for (int dim = 0; dim < KDDIM; dim++) {
			_lowCorner[dim] = low[dim];
			_highCorner[dim] = high[dim];
			_coversWholeDomain[dim] = coversAll[dim];
		}
	}

	/**
	 * compare the tree represented by this node to another tree.
	 */
	bool equals(KDNode& other);

	//! The destructor deletes the childs (recursive call of destructors)
	~KDNode() {
		delete _child1;
		delete _child2;
	}

	/**
	 * @return the area for process rank, i.e. the leaf of this tree with
	 *         (_owningProc == rank) and (_numProcs == 1).
	 *
	 *         If no corresponding node is found, this method returns nullptr!
	 */
	KDNode* findAreaForProcess(int rank);

	//! @brief create an initial decomposition of the domain represented by this node.
	//!
	//! Build a KDTree representing a simple initial domain decomposition by bipartitioning
	//! the area recursively, always in the dimension with the longest extend.
	void buildKDTree();

	/**
	 * @return true, if the node can be resolved for its number of processes (_numProcs),
	 *         i.e. each process can have a subdomain of at least 2 cells per dimension.
	 */
	bool isResolvable();

	/**
	 * @return maximum number of processes, which could be assigned to this node.
	 */
	unsigned int getNumMaxProcs();

	double calculateAvgLoadPerProc() const {
		return _load / ((double) _numProcs);
	}

	void calculateExpectedDeviation(std::vector<double>* accumulatedProcessorSpeeds = nullptr) {
		double meanProcessorSpeed[] = { 1., 1. };
		double averagedMeanProcessorSpeed = 1.;
		if (accumulatedProcessorSpeeds != nullptr && accumulatedProcessorSpeeds->size() != 0) {
			meanProcessorSpeed[0] = ((*accumulatedProcessorSpeeds)[_child2->_owningProc]
					- (*accumulatedProcessorSpeeds)[_owningProc]) / (_child1->_numProcs);
			meanProcessorSpeed[1] = ((*accumulatedProcessorSpeeds)[_child2->_owningProc + _child2->_numProcs]
					- (*accumulatedProcessorSpeeds)[_child2->_owningProc]) / (_child2->_numProcs);
			averagedMeanProcessorSpeed = (meanProcessorSpeed[0] + meanProcessorSpeed[1]) / 2;
		}
		double child1Dev = _child1->calculateAvgLoadPerProc()
				- _optimalLoadPerProcess * meanProcessorSpeed[0] / averagedMeanProcessorSpeed;
		child1Dev = child1Dev * child1Dev;
		double child2Dev = _child2->calculateAvgLoadPerProc()
				- _optimalLoadPerProcess * meanProcessorSpeed[1] / averagedMeanProcessorSpeed;
		child2Dev = child2Dev * child2Dev;
		_expectedDeviation = child1Dev * (double) _child1->_numProcs + child2Dev * (double) _child2->_numProcs;
	}

	void calculateDeviation(std::vector<double>* processorSpeeds = nullptr, const double &totalMeanProcessorSpeed = 1.) {
		if (_numProcs == 1) {
			double speed = 1.;
			if (processorSpeeds != nullptr && processorSpeeds->size() > (unsigned int)_owningProc) {
				speed = (*processorSpeeds)[_owningProc];
			}
			//_deviation = _load - _optimalLoadPerProcess;
			double dev = _load - _optimalLoadPerProcess * speed / totalMeanProcessorSpeed;
			_deviation = dev * dev;
		} else {
			_deviation = _child1->_deviation + _child2->_deviation;
		}
	}

	/**
	 *
	 * @param low
	 * @param high
	 * @param procIDs
	 * @param neighbHaloAreas
	 */
	void getOwningProcs(const int low[KDDIM], const int high[KDDIM], std::vector<int>& procIDs, std::vector<int>& neighbHaloAreas) const;


	/**
	 * Split this node, i.e. create two children (note, that its children must be
	 * nullptr before this call!).
	 *
	 * @param dimension the dimension \in [0;KDDIM-1] along which this node is split
	 * @param splitIndex the index of the corner cell for the new left child
	 *        (note: must be in ] _lowCorner[dimension]; _highCorner[dimension] [.
	 * @param numProcsLeft the number of processors for the left child. The number
	 *        of processors for the right child is calculated.
	 */
	void split(int divDimension, int splitIndex, int numProcsLeft);

	//! @brief prints this (sub-) tree to stdout
	//!
	//! For each node, it is printed whether it is a "LEAF" or a "INNER" node,
	//! The order of printing is a depth-first walk through the tree, children
	//! are always indented two spaces more than there parents
	//! @param prefix A string which is printed in front of each line
	void printTree(const std::string& prefix = "");

	/**
	 * Write the tree represented by this (root-)node to a (binary) file, in order
	 * to be able to restore the decomposition from disk.
	 */
	void serialize(const std::string& fileName);

	/**
	 * Read the tree represented by this (root-)node from a (binary) file.
	 */
	void deserialize(const std::string& fileName);

	/**
	 * plot the leafs of the KDTree with vtk.
	 */
	void plotNode(const std::string& vtkFile, const std::vector<double>* processorSpeeds=nullptr) const;

	/**
	 * Initialize the mpi datatype. Has to be called once initially.
	 */
	static void initMPIDataType() {
		MPIKDNode::initDatatype();
	}

	/**
	 * Free the mpi datatype
	 */
	static void shutdownMPIDataType() {
		MPIKDNode::shutdownDatatype();
	}

	/**
	 * Get a MPIDKNode object representing this KDNode.
	 */
	MPIKDNode getMPIKDNode();

	//! number of procs which share this area
	int _numProcs;
	//! in cells relative to global domain
	int _lowCorner[3];
	//! in cells relative to global domain
	int _highCorner[3];
	//! true if the domain in the given dimension is not divided into more than one process
	bool _coversWholeDomain[KDDIM];

	//! ID of this KDNode 
	int _nodeID;
	//! process which owns this KDNode (only possible for leaf nodes) 
	int _owningProc; // only used if the node is a leaf

	//! "left" child of this KDNode (only used if the child is no leaf)
	KDNode* _child1;
	//! "left" child of this KDNode (only used if the child is no leaf)
	KDNode* _child2;

	double _load;
	double _optimalLoadPerProcess;
	double _expectedDeviation;
	double _deviation;

	// level of this node (at root node, level = 0)
	int _level;

private:

	void serialize(std::ostream& file);

	void deserialize(std::istream& file);

	void plotNode(VTKGridWriterImplementation& writer) const;

};

#endif /* KDNODE_H_ */
