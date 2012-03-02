#ifndef KDNODE_H_
#define KDNODE_H_

// include because of macro KDDIM
#include "parallel/KDDecomposition.h"
/*
 * #include <assert.h>

#define assertion(x) assert(x)
#define assertion1(x,y) assert(x)
 */
#include "parallel/MPIKDNode.h"

//! @brief represents a node in the decomposition tree when using KDDecomposition
//! @author Martin Buchholz
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

	KDNode() : _child1(NULL), _child2(NULL) {
	}

	KDNode(int numP, const int low[KDDIM], const int high[KDDIM], int id, int owner, bool coversAll[KDDIM])
	: _numProcs(numP), _nodeID(id), _owningProc(owner),
	  _child1(NULL), _child2(NULL), _load(0.0), _optimalLoadPerProcess(0.0) {
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
	 *         If no corresponding node is found, this method returns NULL!
	 */
	KDNode* findAreaForProcess(int rank);

	//! @brief create an initial decomposition of the domain represented by this node.
	//!
	//! Build a KDTree representing a simple initial domain decomposition by bipartitioning
	//! the area recursively, always in the dimension with the longest extend.
	void buildKDTree();


	//! @brief prints this (sub-) tree to stdout
	//!
	//! For each node, it is printed whether it is a "LEAF" or a "INNER" node,
	//! The order of printing is a depth-first walk through the tree, children
	//! are always indented two spaces more than there parents
	//! @param prefix A string which is printed in front of each line
	void printTree(std::string prefix = "");

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
};

#endif /*KDNODE_H_*/
