/*
 * KDNode.cpp
 *
 *  Created on: Feb 24, 2012
 *      Author: eckhardw
 */
#include "KDNode.h"

KDNode* KDNode::findAreaForProcess(int rank) {
	if (_numProcs == 1) {
		if (rank == _owningProc) {
			return this;
		} else {
			return NULL;
		}
	}

	if (rank < _child1->_owningProc + _child1->_numProcs) {
		return _child1->findAreaForProcess(rank);
	} else {
		return _child2->findAreaForProcess(rank);
	}
}

bool KDNode::equals(KDNode& other) {
	bool equal = true;

	equal = equal && (_numProcs == other._numProcs);
	equal = equal && (_nodeID == other._nodeID);
	equal = equal && (_owningProc == other._owningProc);

	for (int i = 0; i < KDDIM; i++) {
		equal = equal && (_lowCorner[i] == other._lowCorner[i]);
		equal = equal && (_highCorner[i] == other._highCorner[i]);
		equal = equal && (_coversWholeDomain[i] == other._coversWholeDomain[i]);
	}

	if (_child1 != NULL && other._child1 != NULL) {
		equal == equal && _child1->equals(*other._child1);
	} else {
		equal == equal && (_child1 == other._child1);
	}

	if (_child2 != NULL && other._child2 != NULL) {
		equal == equal && _child2->equals(*other._child2);
	} else {
		equal == equal && (_child2 == other._child2);
	}

	return equal;
}


void KDNode::buildKDTree() {

	if (_numProcs == 1) {
		return;
	}

	bool coversAll[KDDIM];
	int cellsPerDim[KDDIM];
	for (int dim = 0; dim < KDDIM; dim++) {
		coversAll[dim] = _coversWholeDomain[dim];
		cellsPerDim[dim] = _highCorner[dim] - _lowCorner[dim] + 1;
	}
	int divDir = 0;
	int maxCells = cellsPerDim[0];
	if (cellsPerDim[1] > maxCells) {
		divDir = 1;
		maxCells = cellsPerDim[1];
	}
	if (cellsPerDim[2] > maxCells) {
		divDir = 2;
		maxCells = cellsPerDim[2];
	}

	coversAll[divDir] = false;

	int low1[KDDIM];
	int low2[KDDIM];
	int high1[KDDIM];
	int high2[KDDIM];
	int numProcs1;
	int numProcs2;
	int id1;
	int id2;
	int owner1;
	int owner2;

	for (int dim = 0; dim < KDDIM; dim++) {
		low1[dim] = _lowCorner[dim];
		low2[dim] = _lowCorner[dim];
		high1[dim] = _highCorner[dim];
		high2[dim] = _highCorner[dim];
	}

	low1[divDir] = _lowCorner[divDir];
	high1[divDir] = (_highCorner[divDir] + _lowCorner[divDir]) / 2;
	low2[divDir] = high1[divDir] + 1;
	high2[divDir] = _highCorner[divDir];

	numProcs1 = _numProcs / 2;
	numProcs2 = _numProcs - numProcs1;
	id1 = _nodeID + 1;
	id2 = _nodeID + 2 * numProcs1;
	owner1 = _owningProc;
	owner2 = owner1 + numProcs1;
	_child1 = new KDNode(numProcs1, low1, high1, id1, owner1, coversAll);
	_child2 = new KDNode(numProcs2, low2, high2, id2, owner2, coversAll);

	if (numProcs1 > 1) {
		_child1->buildKDTree();
	}

	if (numProcs2 > 1) {
		_child2->buildKDTree();
	}
}
