#include "KDNode.h"

#include <bitset>

#ifdef VTK
#include "io/vtk/VTKGridWriterImplementation.h"
#include "io/vtk/VTKGridVertex.h"
#include "io/vtk/VTKGridCell.h"
#endif
#include "utils/Logger.h"

using namespace Log;

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
	equal = equal && (_level == other._level);

	for (int i = 0; i < KDDIM; i++) {
		equal = equal && (_lowCorner[i] == other._lowCorner[i]);
		equal = equal && (_highCorner[i] == other._highCorner[i]);
		equal = equal && (_coversWholeDomain[i] == other._coversWholeDomain[i]);
	}

	if (_child1 != NULL && other._child1 != NULL) {
		bool childEqual = _child1->equals(*other._child1);
		equal = equal && childEqual;
	} else {
		bool childEqual = (_child1 == other._child1);
		equal = equal && childEqual;
	}

	if (_child2 != NULL && other._child2 != NULL) {
		bool childEqual = _child2->equals(*other._child2);
		equal = equal && childEqual;
	} else {
		bool childEqual = (_child2 == other._child2);
		equal = equal && childEqual;
	}

	return equal;
}


void KDNode::buildKDTree() {

	if (_numProcs == 1) {
		return;
	}

	int cellsPerDim[KDDIM];
	for (int dim = 0; dim < KDDIM; dim++) {
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

	int divIndex = (_highCorner[divDir] + _lowCorner[divDir]) / 2;
	split(divDir, divIndex, _numProcs / 2);

	if (_child1->_numProcs > 1) {
		_child1->buildKDTree();
	}

	if (_child2->_numProcs > 1) {
		_child2->buildKDTree();
	}
}


void KDNode::printTree(std::string prefix) {
// use std::cout as I want to have all nodes at all processes printed
	if (_numProcs == 1) {
		std::cout << prefix << "LEAF: " << _nodeID << ", Owner: " << _owningProc
				<< ", Corners: (" << _lowCorner[0] << ", " << _lowCorner[1] << ", " << _lowCorner[2] << ") / ("
				<< _highCorner[0] << ", " << _highCorner[1] << ", " << _highCorner[2] << "), Load: " << _load << ", OptLoad:" << _optimalLoadPerProcess << " level=" << _level << std::endl;
	}
	else {
		std::cout << prefix << "INNER: " << _nodeID << ", Owner: " << _owningProc
				<< "(" << _numProcs << " procs)" << ", Corners: (" << _lowCorner[0]
				<< ", " << _lowCorner[1] << ", " << _lowCorner[2] << ") / (" << _highCorner[0]
				<< ", " << _highCorner[1] << ", " << _highCorner[2] << ")"
				" child1: " << _child1 << " child2: " << _child2 << ", Load: " << _load << ", OptLoad:" << _optimalLoadPerProcess << " level=" << _level << std::endl;
		std::stringstream childprefix;
		childprefix << prefix << "  ";
		if (_child1 != NULL) {
			_child1->printTree(childprefix.str());
		}
		if (_child2 != NULL) {
			_child2->printTree(childprefix.str());
		}
	}
}

KDNode::MPIKDNode KDNode::getMPIKDNode() {
	std::bitset<3> coversWholeDomain;
	for (int i = 0; i < KDDIM; i++) {
		coversWholeDomain.set(i, _coversWholeDomain[i]);
	}

	int nextSender = _owningProc + 1;
	int leftChildID = -1;
	int rightChildID = -1;

	if (_numProcs > 1) {
		nextSender = _owningProc;
		leftChildID = _child1->_nodeID;
		rightChildID = _child2->_nodeID;
	}

	// const std::bitset<3>& coversWholeDomain, const int& numProcs,
	// const std::vector<int>& lowCorner, const std::vector<int>& highCorner,
	// const int& nodeID, const int& owningProc, const int& firstChildID, const int& secondChildID,
	// const int& nextSendingProcess, const double& load, const double& OptimalLoadPerProcess
	return MPIKDNode(coversWholeDomain, _numProcs, _lowCorner, _highCorner,
				_nodeID, _owningProc, leftChildID, rightChildID, nextSender, _load, _optimalLoadPerProcess,
				_expectedDeviation, _deviation, _level);
}

bool KDNode::isResolvable() {
	int maxProcs = getNumMaxProcs();
	if (maxProcs < _numProcs) {
		return false;
	} else {
		return true;
	}
}

unsigned int KDNode::getNumMaxProcs() {
	unsigned int maxProcs = 1;
	for (int dim = 0; dim < 3; dim++) {
		maxProcs *= (_highCorner[dim] - _lowCorner[dim] + 1) / 2;
	}
	return maxProcs;
}

void KDNode::split(int divDimension, int splitIndex, int numProcsLeft) {
	assert(_numProcs > 1);
	assert(splitIndex > _lowCorner[divDimension]);
	assert(splitIndex < _highCorner[divDimension]);

	bool coversAll[KDDIM];
	for (int dim = 0; dim < KDDIM; dim++) {
		coversAll[dim] = _coversWholeDomain[dim];
	}

	coversAll[divDimension] = false;

	int low1[KDDIM];
	int low2[KDDIM];
	int high1[KDDIM];
	int high2[KDDIM];
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

	low1[divDimension] = _lowCorner[divDimension];
	high1[divDimension] = splitIndex;
	low2[divDimension] = splitIndex + 1;
	high2[divDimension] = _highCorner[divDimension];

	id1 = _nodeID + 1;
	id2 = _nodeID + 2 * numProcsLeft;
	owner1 = _owningProc;
	owner2 = owner1 + numProcsLeft;
	_child1 = new KDNode(numProcsLeft, low1, high1, id1, owner1, coversAll, _level+1);
	_child1->_optimalLoadPerProcess = _optimalLoadPerProcess;
	_child2 = new KDNode(_numProcs - numProcsLeft, low2, high2, id2, owner2, coversAll, _level+1);
	_child2->_optimalLoadPerProcess = _optimalLoadPerProcess;
}

//#define BINARY

void KDNode::serialize(const std::string& fileName) {
	std::ofstream out;
#ifdef BINARY
	out.open(fileName.c_str(), std::ios::binary);
#else
	out.open(fileName.c_str());
#endif
	serialize(out);
	out.close();
}


void KDNode::serialize(std::ostream& file) {
#ifdef BINARY
	file.write((char*) &_numProcs, sizeof(_numProcs));
	file.write((char*) _lowCorner, sizeof(_lowCorner));
	file.write((char*) _highCorner, sizeof(_highCorner));
	file.write((char*) _coversWholeDomain, sizeof(_coversWholeDomain));
	file.write((char*) &_nodeID, sizeof(_nodeID));
	file.write((char*) &_owningProc, sizeof(_owningProc));
	file.write((char*) &_level, sizeof(_level));
#else
	file << _numProcs << " ";
	file << _lowCorner[0] << " " << _lowCorner[1] << " " << _lowCorner[2] << " ";
	file << _highCorner[0] << " " << _highCorner[1] << " " << _highCorner[2] << " ";
	file << _coversWholeDomain[0] << " " << _coversWholeDomain[1] << " " << _coversWholeDomain[2] << " ";
	file << _nodeID << " ";
	file << _owningProc << " ";
	file << _level << " ";
#endif

	if (_numProcs > 1) {
		_child1->serialize(file);
		_child2->serialize(file);
	}
}

void KDNode::deserialize(const std::string& fileName) {
	std::ifstream in;
#ifdef BINARY
	in.open(fileName.c_str(), std::ios::binary);
#else
	in.open(fileName.c_str());
#endif
	deserialize(in);
	in.close();
}

void KDNode::deserialize(std::istream& file) {
	assert(!file.eof());
#ifdef BINARY
	file.read((char*) &_numProcs, sizeof(_numProcs));
	file.read((char*) _lowCorner, sizeof(_lowCorner));
	file.read((char*) _highCorner, sizeof(_highCorner));
	file.read((char*) _coversWholeDomain, sizeof(_coversWholeDomain));
	file.read((char*) &_nodeID, sizeof(_nodeID));
	file.read((char*) &_owningProc, sizeof(_owningProc));
	file.read((char*) &_level, sizeof(_level));
# else
	file >>  _numProcs;
	file >> _lowCorner[0] >> _lowCorner[1] >> _lowCorner[2];
	file >> _highCorner[0] >> _highCorner[1] >> _highCorner[2];
	file >> _coversWholeDomain[0] >> _coversWholeDomain[1] >> _coversWholeDomain[2];
	file >> _nodeID;
	file >> _owningProc;
	file >> _level;
#endif
	if (_numProcs > 1) {
		_child1 = new KDNode();
		_child2 = new KDNode();
		_child1->deserialize(file);
		_child2->deserialize(file);
	}
}

void KDNode::plotNode(const std::string& vtkFile) const {
#ifdef VTK
	VTKGridWriterImplementation writer(_owningProc);
	writer.initializeVTKFile();
	plotNode(writer);
	writer.writeVTKFile(vtkFile);
#else
	global_log->warning() << "KDNode::plotNode() requires vtk output. Compile with -DVTK!"<< std::endl;
#endif
}

void KDNode::plotNode(VTKGridWriterImplementation& writer) const {
#ifdef VTK
	if (_numProcs > 1) {
		_child1->plotNode(writer);
		_child2->plotNode(writer);
	} else {
		VTKGridVertex vertices[8];
		vertices[0].setCoordinates(_lowCorner[0], _lowCorner[1], _lowCorner[2]);
		vertices[1].setCoordinates(_highCorner[0], _lowCorner[1], _lowCorner[2]);
		vertices[2].setCoordinates(_lowCorner[0], _highCorner[1], _lowCorner[2]);
		vertices[3].setCoordinates(_highCorner[0], _highCorner[1], _lowCorner[2]);
		vertices[4].setCoordinates(_lowCorner[0], _lowCorner[1], _highCorner[2]);
		vertices[5].setCoordinates(_highCorner[0], _lowCorner[1], _highCorner[2]);
		vertices[6].setCoordinates(_lowCorner[0], _highCorner[1], _highCorner[2]);
		vertices[7].setCoordinates(_highCorner[0], _highCorner[1], _highCorner[2]);

		VTKGridCell cell;
		cell.setCellData(_owningProc, _load, _level);

		for (int i = 0; i < 8; i++) {
			cell.setVertex(i, &vertices[i]);
		}
		writer.plotCell(cell);
	}
#else
	global_log->warning() << "KDNode::plotNode() requires vtk output. Compile with -DVTK!" << std::endl;
#endif
}
