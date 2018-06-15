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
	if (!currentCell.isHaloCell()) {
		cellProcessor.processCell(currentCell);
		// loop over all forward neighbours
		for (auto& neighbourOffset : this->_forwardNeighbourOffsets) {
			CellTemplate& neighbourCell = this->_cells->at(cellIndex + neighbourOffset);

			const bool sumAllMacroscopic = true;
			cellProcessor.processCellPair(currentCell, neighbourCell, sumAllMacroscopic);
		}
	}
}

#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_HALFSHELLTRAVERSAL_H_ */
