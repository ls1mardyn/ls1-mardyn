/*
 * LinkedCellsHS.cpp
 *
 *  Created on: 27.04.2017
 *      Author: sauermann
 */

#include "particleContainer/LinkedCellsHS.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "utils/mardyn_assert.h"

// Traverse single cell (will be called when not parallelized with OpenMP)
void LinkedCellsHS::traverseCell(long int cellIndex, CellProcessor& cellProcessor){

	ParticleCell& currentCell = _cells[cellIndex];
	if (currentCell.isInnerCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto neighbourOffsetsIter = _forwardNeighbourOffsets.begin();
				neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
			ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
	}
	// Calculating forces between halo cells is unnecessary?
	/*else if (currentCell.isHaloCell()) {
		cellProcessor.processCell(currentCell);
		for (auto neighbourOffsetsIter = _forwardNeighbourOffsets.begin();
				neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
			long int neighbourCellIndex = cellIndex + *neighbourOffsetsIter;
			if ((neighbourCellIndex < 0) || (neighbourCellIndex >= (int) ((_cells.size()))))
				continue;

			ParticleCell& neighbourCell = _cells[neighbourCellIndex];
			if (!neighbourCell.isHaloCell())
				continue;

			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
	} */
	else
	// loop over all boundary cells and calculate forces to forward and backward neighbours
	if (currentCell.isBoundaryCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto neighbourOffsetsIter = _forwardNeighbourOffsets.begin();
				neighbourOffsetsIter != _forwardNeighbourOffsets.end(); neighbourOffsetsIter++) {
			ParticleCell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
	} // if ( isBoundaryCell() )

}

// Traverse all cells with OpenMP
void LinkedCellsHS::traverseCellsC08(CellProcessor& cellProcessor) {
	//TODO ___Implement me
	std::cout << "Not yet implemented\n\n\n\n\n";
 	mardyn_exit(1);
}
