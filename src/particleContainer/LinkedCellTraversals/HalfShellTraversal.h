/*
 * HalfShellTraversal.h
 *
 *  Created on: 08.06.2017
 *      Author: sauermann
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_HALFSHELLTRAVERSAL_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_HALFSHELLTRAVERSAL_H_

#include "particleContainer/LinkedCellTraversals/OriginalCellPairTraversal.h"

struct HalfShellTraversalData : OriginalCellPairTraversalData {
};

template <class CellTemplate>
class HalfShellTraversal: public OriginalCellPairTraversal<CellTemplate> {
	using OriginalCellPairTraversal<CellTemplate>::OriginalCellPairTraversal; // Inheriting constructors

	bool requiresForceExchange() const override {return true;}
public:

	void processBaseCell(CellProcessor& cellProcessor, unsigned long cellIndex) const override;
};

template<class CellTemplate>
inline void HalfShellTraversal<CellTemplate>::processBaseCell(CellProcessor& cellProcessor, unsigned long cellIndex) const{

	CellTemplate& currentCell = this->_cells->at(cellIndex);
	if (currentCell.isInnerCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto& neighbourOffset : this->_forwardNeighbourOffsets) {
			CellTemplate& neighbourCell = this->_cells->at(cellIndex + neighbourOffset);
			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
	}
	else
	// loop over all boundary cells and calculate forces to forward neighbours
	if (currentCell.isBoundaryCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto& neighbourOffset : this->_forwardNeighbourOffsets) {
			CellTemplate& neighbourCell = this->_cells->at(cellIndex + neighbourOffset);
			cellProcessor.processCellPair(currentCell, neighbourCell);
		}
	} // if ( isBoundaryCell() )

}

#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_HALFSHELLTRAVERSAL_H_ */
