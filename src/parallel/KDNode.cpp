#include "KDNode.h"

#include <algorithm>
#include <bitset>

#ifdef VTK
#include "io/vtk/VTKGridWriterImplementation.h"
#include "io/vtk/VTKGridVertex.h"
#include "io/vtk/VTKGridCell.h"
#endif
#include "KDDStaticValues.h"
#include "utils/Logger.h"

#include <algorithm> /* for min and max ?*/

using namespace Log;

KDNode* KDNode::findAreaForProcess(int rank) {
	if (_numProcs == 1) {
		if (rank == _owningProc) {
			return this;
		} else {
			return nullptr;
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

	if (_child1 && other._child1 ) {
		bool childEqual = _child1->equals(*other._child1);
		equal = equal && childEqual;
	} else {
		bool childEqual = (_child1 == other._child1);
		equal = equal && childEqual;
	}

	if (_child2 && other._child2) {
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

	// split in the dimension with the most cells
	std::vector<int> cellsPerDim(KDDIM);
	for (int dim = 0; dim < KDDIM; dim++) {
		cellsPerDim[dim] = _highCorner[dim] - _lowCorner[dim] + 1;
	}
	auto max_element_ptr = std::max_element(cellsPerDim.begin(), cellsPerDim.end());
	int divDir = max_element_ptr - cellsPerDim.begin();
	int divIndex = (_highCorner[divDir] + _lowCorner[divDir]) / 2;
	split(divDir, divIndex, _numProcs / 2);

	if (_child1->_numProcs > 1) {
		_child1->buildKDTree();
	}
	if (_child2->_numProcs > 1) {
		_child2->buildKDTree();
	}
}

void KDNode::printTree(const std::string& prefix) {
	// use std::cout as I want to have all nodes at all processes printed
	if (_numProcs == 1) {
		std::cout << prefix << "LEAF: " << _nodeID << ", Owner: " << _owningProc << ", Corners: (" << _lowCorner[0]
				  << ", " << _lowCorner[1] << ", " << _lowCorner[2] << ") / (" << _highCorner[0] << ", "
				  << _highCorner[1] << ", " << _highCorner[2]
				  << "), Length(cells): " << _highCorner[0] - _lowCorner[0] + 1 << "x"
				  << _highCorner[1] - _lowCorner[1] + 1 << "x" << _highCorner[2] - _lowCorner[2] + 1
				  << " Load: " << _load << ", deviation: " << _deviation << ", OptLoad:" << _optimalLoadPerProcess
				  << " level=" << _level << std::endl;
	} else {
		std::cout << prefix << "INNER: " << _nodeID << ", Owner: " << _owningProc << "(" << _numProcs << " procs)"
				  << ", Corners: (" << _lowCorner[0] << ", " << _lowCorner[1] << ", " << _lowCorner[2] << ") / ("
				  << _highCorner[0] << ", " << _highCorner[1] << ", " << _highCorner[2] << ") child1: " << _child1
				  << " child2: " << _child2 << ", Load: " << _load << ", deviation: " << _deviation
				  << ", OptLoad:" << _optimalLoadPerProcess << " level=" << _level << std::endl;
		std::stringstream childprefix;
		childprefix << prefix << "  ";
		if (_child1) {
			_child1->printTree(childprefix.str());
		}
		if (_child2) {
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
				_nodeID, _owningProc, leftChildID, rightChildID, nextSender, _load, _optimalLoadPerProcess, _deviationLowerBound, _deviation, _level);
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
	// we need at least KDDStaticValues::minNumCellsPerDimension cells in each dimension per process.
	unsigned int maxProcs = 1;
	for (int dim = 0; dim < 3; dim++) {
		maxProcs *= (_highCorner[dim] - _lowCorner[dim] + 1) / KDDStaticValues::minNumCellsPerDimension;
	}
	return maxProcs;
}

void KDNode::split(int divDimension, int splitIndex, int numProcsLeft) {
	mardyn_assert(_numProcs > 1);
	mardyn_assert(splitIndex >= _lowCorner[divDimension] + (KDDStaticValues::minNumCellsPerDimension - 1));
	mardyn_assert(splitIndex < _highCorner[divDimension]);

	bool coversAll[KDDIM];
	for (int dim = 0; dim < KDDIM; dim++) {
		coversAll[dim] = _coversWholeDomain[dim];
	}

	// if we cut it along the dimension divDimension, it can no longer cover the whole domain
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
	mardyn_assert(!file.eof());
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

void KDNode::plotNode(const std::string& vtkFile, const std::vector<double>* processorSpeeds) const {
#ifdef VTK
	VTKGridWriterImplementation writer(_owningProc, processorSpeeds);
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
		cell.setIndex(_nodeID);
		for (int i = 0; i < 8; i++) {
			cell.setVertex(i, &vertices[i]);
		}
		writer.plotCell(cell);
	}
#else
	global_log->warning() << "KDNode::plotNode() requires vtk output. Compile with -DVTK!" << std::endl;
#endif
}

void KDNode::getOwningProcs(const int low[KDDIM], const int high[KDDIM], std::vector<int>& procIDs, std::vector<int>& neighbAreas) const {

	int overlapLow[3];
	int overlapHigh[3];

	bool areasIntersect = true;
	for (int dim = 0; dim < KDDIM; dim++) {
		overlapLow[dim]  = std::max(low[dim],  _lowCorner[dim]);
		overlapHigh[dim] = std::min(high[dim], _highCorner[dim]);
		areasIntersect &= overlapLow[dim] <= _highCorner[dim] and overlapHigh[dim] >= _lowCorner[dim];
	}

	if (areasIntersect) {
		if (_numProcs > 1) {
			_child1->getOwningProcs(overlapLow, overlapHigh, procIDs, neighbAreas);
			_child2->getOwningProcs(overlapLow, overlapHigh, procIDs, neighbAreas);
		}
		else { // leaf node found, add it (and it's area)
			procIDs.push_back(_owningProc);
			for (int d = 0; d < KDDIM; ++d) {
				neighbAreas.push_back(overlapLow[d]);
				neighbAreas.push_back(overlapHigh[d]);
			}
		}
	}
}
