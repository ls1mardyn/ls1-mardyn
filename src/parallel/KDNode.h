#ifndef KDNODE_H_
#define KDNODE_H_

//! @brief represents a node in the decomposition tree when using KDDecomposition
//! @author Martin Buchholz
//! 
//! The KDDecomposition decomposes the domain by recursively splitting the domain
//! into smaller parts. This class is used to represent this decomposition.
//! The root node of the decomposition covers the whole domain, in the first
//! splitting step, this domain is divided into two parts, where each part
//! then has to be divided into several smaller parts. How many parts/regions there
//! will be depends on the number of processes, each process will get one region.
//! So the KNNode also has to store how many process "share" the current region
//! The leaf nodes of the tree represent the region of the single processes.
//! The regions (the size of the regions) is not stored in some floating-point length unit
//! but in cells. So it is assumed that the domain is discretised with cells and
//! the decomposition is based on distributing those cells (blocks of cells) to the
//! processes
class KDNode {
public:
	//! This constructor just sets the child pointers to NULL
	KDNode() {
		_child1 = NULL;
		_child2 = NULL;
	}

	//! This constructor sets all member variables to the given values (child pointers to NULL)
	KDNode(int numP, int low[KDDIM], int high[KDDIM], int id, int owner, bool coversAll[KDDIM]) {
		_numProcs = numP;
		for (int dim = 0; dim < KDDIM; dim++) {
			_lowCorner[dim] = low[dim];
			_highCorner[dim] = high[dim];
			_coversWholeDomain[dim] = coversAll[dim];
		}
		_nodeID = id;
		_owningProc = owner;
		_child1 = NULL;
		_child2 = NULL;
	}

	//! The destructor deletes the childs (recursive call of destructors)
	~KDNode() {
		delete _child1;
		delete _child2;
	}

	//! checks whether the given cell position belongs to the region of this KDNode
	bool cellBelongsToRegion(int cellPos[KDDIM]) {
		bool isInside = true;
		for (int dim = 0; dim < KDDIM; dim++) {
			if (cellPos[dim] < _lowCorner[dim] || cellPos[dim] > _highCorner[dim])
				isInside = false;
		}
		return isInside;
	}

	//! checks whether the given cell position is in the boundary region of this KDNode 
	bool regionCellBelongsToBoundary(int cellPos[KDDIM]) {
		bool isBoundary = false;
		for (int dim = 0; dim < KDDIM; dim++) {
			if (cellPos[dim] == _lowCorner[dim] || cellPos[dim] == _highCorner[dim])
				isBoundary = true;
		}
		return isBoundary;
	}

	//! number of procs which share this area
	int _numProcs;
	//! in cells relative to global domain
	int _lowCorner[KDDIM];
	//! in cells relative to global domain
	int _highCorner[KDDIM];
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
};

#endif /*KDNODE_H_*/
