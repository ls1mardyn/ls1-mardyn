/*
 * LeafNodesContainer.cpp
 *
 *  Created on: Feb 6, 2015
 *      Author: tchipev
 */

#include "LeafNodesContainer.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "bhfmm/cellProcessors/VectorizedChargeP2PCellProcessor.h"
#include "bhfmm/cellProcessors/SimpleCellProcessor.h"
#include "WrapOpenMP.h"

#include <cmath>
#include <bhfmm/FastMultipoleMethod.h>

using namespace std;
using Log::global_log;

namespace bhfmm {

LeafNodesContainer::LeafNodesContainer(double bBoxMin[3],
									   double bBoxMax[3],
									   double LJCellLength[3],
									   unsigned subdivisionFactor,
									   bool periodic
#ifdef QUICKSCHED
									   , qsched* scheduler ) : _scheduler(scheduler) {
#else
									   ) {
#endif

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
	calculateCellPairOffsets();

	// get number of active threads
	int num_active_threads = mardyn_get_max_threads();

	int strides[3];
	if(num_active_threads > 1) {
		strides[0] = 2;
		strides[1] = 2;
		strides[2] = 2;
		_numActiveColours = 8;
	} else {
		strides[0] = 1;
		strides[1] = 1;
		strides[2] = 1;
		_numActiveColours = 1;
	}
	_cellIndicesPerColour.resize(_numActiveColours);

	// initialize the vector of cells per color
	// with the index of each cell in the respective color
	// in a 8-coloring-scheme each cell in a 2x2x2 block has a unique color
	for (unsigned col = 0; col < _numActiveColours; ++col) {
		int start_indices[3];
		threeDIndexOfCellIndex((int)(col), start_indices, strides);

		_cellIndicesPerColour[col].clear();

		// compute indices first
		for (int z = start_indices[2]; z < _numCellsPerDimension[2]-1; z += strides[2]) {
			for (int y = start_indices[1]; y < _numCellsPerDimension[1]-1; y += strides[1]) {
				for (int x = start_indices[0]; x < _numCellsPerDimension[0]-1; x += strides[0]) {
					long int cellIndex = cellIndexOf3DIndex(x,y,z);
					_cellIndicesPerColour[col].push_back(cellIndex);
				} // x
			} // y
		} // z
	} // col
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
				cellPosition[0] = _boundingBoxMin[0] + (ix  ) * _cellLength[0];
				cellPosition[1] = _boundingBoxMin[1] + (iy  ) * _cellLength[1];
				cellPosition[2] = _boundingBoxMin[2] + (iz  ) * _cellLength[2];
				_cells[cellIndex].setBoxMax(cellPosition);
			}
		}
	}
	generateP2PTasks();
}

void LeafNodesContainer::generateP2PTasks() {
#ifdef QUICKSCHED
	unsigned long cellIndex;
	qsched_res_t  resourceId;
	//TODO: Make taskBlockSize dynamic / autotuning / xml
	struct qsched_payload payload;
	payload.taskBlockSize[0] = 2;
	payload.taskBlockSize[1] = 2;
	payload.taskBlockSize[2] = 2;
	payload.leafNodesContainer = this;


	global_log->info() << "LeafNodesContainer: Generating resource and task ids" << std::endl;
	for (auto z = 0; z < _numCellsPerDimension[2]; ++z) {
		for (auto y = 0; y < _numCellsPerDimension[1]; ++y) {
			for (auto x = 0; x < _numCellsPerDimension[0]; ++x) {
				cellIndex  = cellIndexOf3DIndex(x,y,z);
				resourceId = qsched_addres(_scheduler, qsched_owner_none, qsched_res_none);
				_cells[cellIndex].setResourceId(resourceId);
				// only create tasks with offset blocksize-1.
				// -1 because they need to overlap
				// skip tasks for rear halo layer as they would only contain halo cells
				if ((z % (payload.taskBlockSize[2] - 1) == 0
					 && y % (payload.taskBlockSize[1] - 1) == 0
					 && x % (payload.taskBlockSize[0] - 1) == 0)
					&&
					(x < _numCellsPerDimension[0] - 1
					 && y < _numCellsPerDimension[1] - 1
					 && z < _numCellsPerDimension[2] - 1)) {
					// P2P TASK
					// also save the pointers as long
					payload.cell.coordinates[0] = x;
					payload.cell.coordinates[1] = y;
					payload.cell.coordinates[2] = z;
					_cells[cellIndex].setP2PId(qsched_addtask(_scheduler,
															  FastMultipoleMethod::P2P,
															  task_flag_none,
															  &payload,
															  sizeof(payload),
															  1)
					);
				}
				// Pre and post tasks for every cell
				// PREPROCESS TASK
				payload.cell.pointer = &_cells[cellIndex];
				_cells[cellIndex].setPreprocessId(qsched_addtask(_scheduler,
																 FastMultipoleMethod::PreprocessCell,
																 task_flag_none,
																 &payload,
																 sizeof(payload),
																 1)
				);
				// POSTPROCESS TASK
				_cells[cellIndex].setPostprocessId(qsched_addtask(_scheduler,
																  FastMultipoleMethod::PostprocessCell,
																  task_flag_none,
																  &payload,
																  sizeof(payload),
																  1)
				);
			} /* end for-x */
		} /* end for-y*/
	} /* end for-z */

	// set dependencies
	global_log->info() << "LeafNodesContainer: Setting task dependencies" << std::endl;
	for (auto z = 0; z < _numCellsPerDimension[2] - 1; z += payload.taskBlockSize[2] - 1) {
		for (auto y = 0; y < _numCellsPerDimension[1] - 1; y += payload.taskBlockSize[1] - 1) {
			for (auto x = 0; x < _numCellsPerDimension[0] - 1; x += payload.taskBlockSize[0] - 1) {
				cellIndex = cellIndexOf3DIndex(x, y, z);

				for (auto i = 0; i < payload.taskBlockSize[0]
								 && x + i < _numCellsPerDimension[0]; ++i) {
					for (auto j = 0; j < payload.taskBlockSize[1]
									 && y + j < _numCellsPerDimension[1]; ++j) {
						for (auto k = 0; k < payload.taskBlockSize[2]
										 && z + k < _numCellsPerDimension[2]; ++k) {
							// create locks for only for resources at edges
							if(i == payload.taskBlockSize[0] - 1
							   || j == payload.taskBlockSize[1] - 1
							   || k == payload.taskBlockSize[2] - 1){
								qsched_addlock(_scheduler,
											   _cells[cellIndex]
													   .getTaskData()
													   ._P2PId,
											   _cells[cellIndexOf3DIndex(x + i, y + j, z + k)]
													   .getTaskData()
													   ._resourceId);
							}
							// one preprocess (partly) unlocks 8 P2P
							qsched_addunlock(_scheduler,
											 _cells[cellIndexOf3DIndex(x + i, y + j, z + k)]
													 .getTaskData()
													 ._preprocessId,
											 _cells[cellIndex]
													 .getTaskData()
													 ._P2PId);
							// every P2P (partly) unlocks one postprocess
							qsched_addunlock(_scheduler,
											 _cells[cellIndex]
													 .getTaskData()
													 ._P2PId,
											_cells[cellIndexOf3DIndex(x + i, y + j, z + k)]
													 .getTaskData()
													 ._postprocessId);
						}
					}
				}
			} /* end for-x */
		} /* end for-y*/
	} /* end for-z */
#endif /* QUICKSCHED */
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
	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (long int cellIndex = 0; cellIndex < (long int) _cells.size(); cellIndex++){
		_cells[cellIndex].removeAllParticles();
	}
}

unsigned long int LeafNodesContainer::getCellIndexOfMolecule(Molecule* molecule) const {
	int cellIndex[3]; // 3D Cell index

	for (int dim = 0; dim < 3; dim++) {
		cellIndex[dim] = (int) floor((molecule->r(dim) - _haloBoundingBoxMin[dim]) / _cellLength[dim]);

	}
	return this->cellIndexOf3DIndex( cellIndex[0], cellIndex[1], cellIndex[2] );
}

void LeafNodesContainer::traverseCells(SimpleCellProcessor& cellProcessor){
	cellProcessor.initTraversal();

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(static)
	#endif
	for (unsigned int cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
		cellProcessor.processCell(_cells[cellIndex]);
	}

	cellProcessor.endTraversal();
}

void LeafNodesContainer::traverseCellPairs(VectorizedChargeP2PCellProcessor& cellProcessor) {
	#ifndef NDEBUG
		global_log->debug() << "LeafNodesContainer::traverseCells: Processing pairs." << endl;
		global_log->debug() << "_minNeighbourOffset=" << _minNeighbourOffset << "; _maxNeighbourOffset=" << _maxNeighbourOffset<< endl;
	#endif

	#if defined(_OPENMP)
		traverseCellPairsC08(cellProcessor);
	#else
		traverseCellPairsOrig(cellProcessor);
	#endif
}

void LeafNodesContainer::traverseCellPairsOrig(VectorizedChargeP2PCellProcessor& cellProcessor) {	vector<unsigned long>::iterator neighbourOffsetsIter;

	cellProcessor.initTraversal();
	// preprocess all cells
	for (unsigned int cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
		cellProcessor.preprocessCell(_cells[cellIndex]);
	}

	// loop over all inner cells and calculate forces to forward neighbours
	for (unsigned int cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
		ParticleCellPointers& currentCell = _cells[cellIndex];

		if (currentCell.isInnerCell()) {
			cellProcessor.processCell(currentCell);
			// loop over all neighbours
			for (neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				ParticleCellPointers& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}
		}

		if (currentCell.isHaloCell()) {
			cellProcessor.processCell(currentCell);
			for (neighbourOffsetsIter = _forwardNeighbourOffsets.begin(); neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				int neighbourCellIndex = cellIndex + *neighbourOffsetsIter;
				if ((neighbourCellIndex < 0) || (neighbourCellIndex >= (int) (_cells.size())))
					continue;
				ParticleCellPointers& neighbourCell = _cells[neighbourCellIndex];
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
				ParticleCellPointers& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
				cellProcessor.processCellPair(currentCell, neighbourCell);
			}

			// loop over all backward neighbours. calculate only forces
			// to neighbour cells in the halo region, all others already have been calculated
			for (neighbourOffsetsIter = _backwardNeighbourOffsets.begin(); neighbourOffsetsIter != _backwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
				ParticleCellPointers& neighbourCell = _cells[cellIndex - *neighbourOffsetsIter];
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

void LeafNodesContainer::traverseCellPairsC08(VectorizedChargeP2PCellProcessor& cellProcessor) {
	#if defined(_OPENMP)
		cellProcessor.initTraversal();

		// preprocess all cells
		#pragma omp parallel for schedule(static)
		for (unsigned int cellIndex = 0; cellIndex < _cells.size(); ++cellIndex) {
			cellProcessor.preprocessCell(_cells[cellIndex]);
		}

		// process all cells
		#pragma omp parallel
		{
			for (unsigned col = 0; col < _numActiveColours; ++col) {
				const int numIndicesOfThisColour = _cellIndicesPerColour[col].size();

				#pragma omp for schedule(dynamic)
				for(int i = 0; i < numIndicesOfThisColour; ++i) {
					long int baseIndex = _cellIndicesPerColour[col][i];

					c08Step(baseIndex, cellProcessor);

				} // for-loop over indices of this colour
			} // for-loop over colours
		} // end pragma omp parallel

		// postprocess all cells
		#pragma omp parallel for schedule(static)
		for (unsigned int cellIndex = 0; cellIndex < _cells.size(); cellIndex++) {
				cellProcessor.postprocessCell(_cells[cellIndex]);
		}

		cellProcessor.endTraversal();
	#endif
}

void LeafNodesContainer::c08Step(long int baseIndex, VectorizedChargeP2PCellProcessor &cellProcessor) {
	const int num_pairs = _cellPairOffsets.size();
	for(int j = 0; j < num_pairs; ++j) {
		pair<long int, long int> current_pair = _cellPairOffsets[j];

		long int offset1 = current_pair.first;
		long int cellIndex1 = baseIndex + offset1;
		if ((cellIndex1 < 0) || (cellIndex1 >= (int) (_cells.size())))
			continue;

		long int offset2 = current_pair.second;
		long int cellIndex2 = baseIndex + offset2;
		if ((cellIndex2 < 0) || (cellIndex2 >= (int) (_cells.size())))
			continue;

		ParticleCellPointers& cell1 = _cells[cellIndex1];
		ParticleCellPointers& cell2 = _cells[cellIndex2];

		if(cell1.isHaloCell() and cell2.isHaloCell()) {
			continue;
		}

		if(cellIndex1 == cellIndex2) {
			cellProcessor.processCell(cell1);
		}
		else {
			if(!cell1.isHaloCell()) {
				cellProcessor.processCellPair(cell1, cell2);
			}
			else {
				cellProcessor.processCellPair(cell2, cell1);
			}
		}
	}
}

LeafNodesContainer::~LeafNodesContainer() {
}

void LeafNodesContainer::calculateCellPairOffsets() {
	_cellPairOffsets.clear();
	_cellPairOffsets.reserve(14);

	long int o   = cellIndexOf3DIndex(0,0,0); // origin
	long int x   = cellIndexOf3DIndex(1,0,0); // displacement to the right
	long int y   = cellIndexOf3DIndex(0,1,0); // displacement ...
	long int z   = cellIndexOf3DIndex(0,0,1);
	long int xy  = cellIndexOf3DIndex(1,1,0);
	long int yz  = cellIndexOf3DIndex(0,1,1);
	long int xz  = cellIndexOf3DIndex(1,0,1);
	long int xyz = cellIndexOf3DIndex(1,1,1);

	// minimize number of cells simultaneously in memory:

	_cellPairOffsets.push_back(make_pair(o, xyz));
	// evict xyz

	_cellPairOffsets.push_back(make_pair(o, yz ));
	_cellPairOffsets.push_back(make_pair(x, yz ));
	// evict yz

	_cellPairOffsets.push_back(make_pair(o, x  ));

	_cellPairOffsets.push_back(make_pair(o, xy ));
	_cellPairOffsets.push_back(make_pair(xy, z ));
	// evict xy

	_cellPairOffsets.push_back(make_pair(o, z  ));
	_cellPairOffsets.push_back(make_pair(x, z  ));
	_cellPairOffsets.push_back(make_pair(y, z  ));
	// evict z

	_cellPairOffsets.push_back(make_pair(o, y  ));
	_cellPairOffsets.push_back(make_pair(x, y  ));
	// evict x

	_cellPairOffsets.push_back(make_pair(o, xz ));
	_cellPairOffsets.push_back(make_pair(y, xz ));
	// evict xz

	_cellPairOffsets.push_back(make_pair(o, o  ));
}

void LeafNodesContainer::threeDIndexOfCellIndex(int ind, int r[3], int dim[3]) const {
	r[2] = ind / (dim[0] * dim[1]);
	r[1] = (ind - r[2] * dim[0] * dim[1]) / dim[0];
	r[0] = ind - dim[0] * (r[1] + dim[1] * r[2]);
}

	const int *LeafNodesContainer::getNumCellsPerDimension() const {
		return _numCellsPerDimension;
	}

	vector<ParticleCellPointers> & LeafNodesContainer::getCells() {
		return _cells;
	}

} /* namespace bhfmm */
