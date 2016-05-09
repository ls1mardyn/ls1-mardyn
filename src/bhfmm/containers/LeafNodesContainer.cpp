/*
 * LeafNodesContainer.cpp
 *
 *  Created on: Feb 6, 2015
 *      Author: tchipev
 */

#include "LeafNodesContainer.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "particleContainer/adapter/CellProcessor.h"


#include <cmath>

using namespace std;
using Log::global_log;

namespace bhfmm {

LeafNodesContainer::LeafNodesContainer(double bBoxMin[3],
		double bBoxMax[3], double LJCellLength[3], unsigned subdivisionFactor,
		bool periodic) {

	_periodicBC = periodic; // a hack-in workaround to disable periodicity is
	// simply to change the addParticle functionality to filter out by
	// the bounding box, instead of by the HaloBoundingBox


	for (int d = 0; d < 3; d++) {
		_boundingBoxMin[d] = bBoxMin[d];
		_boundingBoxMax[d] = bBoxMax[d];
	}

	unsigned totNumCells = 1;

	for (int d = 0; d < 3; ++d) {
		_cellLength[d] = LJCellLength[d] / subdivisionFactor;
		_numInnerCellsPerDimension[d] = floor((_boundingBoxMax[d] - _boundingBoxMin[d]) / _cellLength[d] + 0.5);
		_numCellsPerDimension[d] = _numInnerCellsPerDimension[d] + 2;

		totNumCells *= _numCellsPerDimension[d];

		_haloBoundingBoxMin[d] = _boundingBoxMin[d] - _cellLength[d];
		_haloBoundingBoxMax[d] = _boundingBoxMax[d] + _cellLength[d];

	}

	_cells.reserve(totNumCells);
	_cells.resize(totNumCells);

	initializeCells();
	calculateNeighbourIndices();
}

void LeafNodesContainer::initializeCells() {

	unsigned long cellIndex;
	for (int iz = 0; iz < _numCellsPerDimension[2]; ++iz) {
		for (int iy = 0; iy < _numCellsPerDimension[1]; ++iy) {
			for (int ix = 0; ix < _numCellsPerDimension[0]; ++ix) {

				cellIndex = cellIndexOf3DIndex(ix, iy, iz);
				_cells[cellIndex].skipCellFromHaloRegion();
				_cells[cellIndex].skipCellFromBoundaryRegion();
				_cells[cellIndex].skipCellFromInnerRegion();

				if (ix == 0 || iy == 0 || iz == 0 ||
				    ix == _numCellsPerDimension[0]-1 ||
				    iy == _numCellsPerDimension[1]-1 ||
				    iz == _numCellsPerDimension[2]-1) {
					_cells[cellIndex].assignCellToHaloRegion();
				}
				else if (ix == 1 || iy == 1 || iz == 1 ||
				         ix == _numCellsPerDimension[0]-2 ||
				         iy == _numCellsPerDimension[1]-2 ||
				         iz == _numCellsPerDimension[2]-2) {
					_cells[cellIndex].assignCellToBoundaryRegion();
				}
				else {
					_cells[cellIndex].assignCellToInnerRegion();
				}

				double cellPosition[3];
				cellPosition[0] = _boundingBoxMin[0] + (ix-1) * _cellLength[0];
				cellPosition[1] = _boundingBoxMin[1] + (iy-1) * _cellLength[1];
				cellPosition[2] = _boundingBoxMin[2] + (iz-1) * _cellLength[2];
				_cells[cellIndex].setBoxMin(cellPosition);
			}
		}
	}
}


void LeafNodesContainer::calculateNeighbourIndices() {
	global_log->debug() << "Setting up cell neighbour index lists." << endl;
	_forwardNeighbourOffsets.clear();
	_backwardNeighbourOffsets.clear();
	_maxNeighbourOffset = 0;
	_minNeighbourOffset = 0;
	for (int zIndex = -1; zIndex <= 1; zIndex++) {
		for (int yIndex = -1; yIndex <= 1; yIndex++) {
			for (int xIndex = -1; xIndex <= 1; xIndex++) {
				long int offset = cellIndexOf3DIndex(xIndex, yIndex, zIndex);
				if (offset > 0) {
					_forwardNeighbourOffsets.push_back(abs(offset));
					if (offset > _maxNeighbourOffset) {
						_maxNeighbourOffset = offset;
					}
				}
				if (offset < 0) {
					_backwardNeighbourOffsets.push_back(abs(offset));
					if (abs(offset) > _minNeighbourOffset) {
						_minNeighbourOffset = abs(offset);
					}
				}
			}
		}
	}

	global_log->info() << "Neighbour offsets are bounded by "
			<< _minNeighbourOffset << ", " << _maxNeighbourOffset << endl;
}

long int LeafNodesContainer::cellIndexOf3DIndex(int xIndex, int yIndex, int zIndex) const {
	return (zIndex * _numCellsPerDimension[1] + yIndex) * _numCellsPerDimension[0] + xIndex;
}

void LeafNodesContainer::addParticle(Molecule& particle) {

	bool insert;
	if(_periodicBC == true) {
		insert = particle.inBox(_haloBoundingBoxMin, _haloBoundingBoxMax);
	} else {
		insert = particle.inBox(_boundingBoxMin, _boundingBoxMax);
	}


	if (insert) {
		int cellIndex = getCellIndexOfMolecule(&particle);
		_cells[cellIndex].addParticle(&particle);
	}
}

void LeafNodesContainer::clearParticles() {
	// clear all Cells
	std::vector<ParticleCell>::iterator celliter;
	for (celliter = (_cells).begin(); celliter != (_cells).end(); ++celliter) {
		(*celliter).removeAllParticles();
	}
}

unsigned long int LeafNodesContainer::getCellIndexOfMolecule(Molecule* molecule) const {
	int cellIndex[3]; // 3D Cell index

	for (int dim = 0; dim < 3; dim++) {
		cellIndex[dim] = (int) floor((molecule->r(dim) - _haloBoundingBoxMin[dim]) / _cellLength[dim]);

	}
	return this->cellIndexOf3DIndex( cellIndex[0], cellIndex[1], cellIndex[2] );
}

void LeafNodesContainer::traverseCellPairs(CellProcessor& cellProcessor) {

	vector<unsigned long>::iterator neighbourOffsetsIter;

#ifndef NDEBUG
	global_log->debug() << "LeafNodesContainer::traverseCells: Processing pairs." << endl;
	global_log->debug() << "_minNeighbourOffset=" << _minNeighbourOffset << "; _maxNeighbourOffset=" << _maxNeighbourOffset<< endl;
#endif

	cellProcessor.initTraversal(_maxNeighbourOffset + _minNeighbourOffset +1);
	// preprocess all cells
	for (unsigned int cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
		cellProcessor.preprocessCell(_cells[cellIndex]);
	}

	// loop over all inner cells and calculate forces to forward neighbours
	for (unsigned int cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
		ParticleCell& currentCell = _cells[cellIndex];

		if (currentCell.isInnerCell()) {
			cellProcessor.processCell(currentCell);
			// loop over all neighbours
			for (neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}

		if (currentCell.isHaloCell()) {
			cellProcessor.processCell(currentCell);
			for (neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				int neighbourCellIndex = cellIndex + *neighbourOffsetsIter;
				if ((neighbourCellIndex < 0) || (neighbourCellIndex >= (int) (_cells.size())))
					continue;
				ParticleCell& neighbourCell = _cells[neighbourCellIndex];
				if (!neighbourCell.isHaloCell())
					continue;

				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}

		// loop over all boundary cells and calculate forces to forward and backward neighbours
		if (currentCell.isBoundaryCell()) {
			cellProcessor.processCell(currentCell);

			// loop over all forward neighbours
			for (neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}

			// loop over all backward neighbours. calculate only forces
			// to neighbour cells in the halo region, all others already have been calculated
			for (neighbourOffsetsIter = _backwardNeighbourOffsets.begin(); neighbourOffsetsIter != _backwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				ParticleCell& neighbourCell = _cells[cellIndex - *neighbourOffsetsIter];
				if (neighbourCell.isHaloCell()) {
					cellProcessor.processCellPair(currentCell, neighbourCell);
				}
			}
		} // if ( isBoundaryCell() )
	} // loop over all cells

	// postprocess all cells
	for (unsigned int cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
			cellProcessor.postprocessCell(_cells[cellIndex]);
	}

	cellProcessor.endTraversal();
}

LeafNodesContainer::~LeafNodesContainer() {
}

} /* namespace bhfmm */
