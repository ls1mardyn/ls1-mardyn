/*
 * NeutralTerritoryTraversal.h
 *
 *  Created on: 25.11.2020
 *      Author: seckler
 */

#pragma once

#include <array>
#include "particleContainer/LinkedCellTraversals/CellPairTraversals.h"

struct NeutralTerritoryTraversalData : CellPairTraversalData {};

template <class CellTemplate>
class NeutralTerritoryTraversal : public CellPairTraversals<CellTemplate> {
public:
	NeutralTerritoryTraversal(std::vector<CellTemplate>& cells, const std::array<unsigned long, 3>& dims, double cellLength[3], double cutoff)
		: CellPairTraversals<CellTemplate>(cells, dims) {
        computeOffsets3D(cellLength, cutoff);
        computeOffsets();
	}
	~NeutralTerritoryTraversal() override = default;

	/**
	 * Reset all necessary data.
	 */
	void rebuild(std::vector<CellTemplate>& cells, const std::array<unsigned long, 3>& dims, double cellLength[3], double cutoff,
				 CellPairTraversalData* data) override {
		CellPairTraversals<CellTemplate>::rebuild(cells, dims, cellLength, cutoff, data);
        computeOffsets3D(cellLength, cutoff);
        computeOffsets();
	}

	void traverseCellPairs(CellProcessor& cellProcessor) override;
	void traverseCellPairsOuter(CellProcessor& cellProcessor) override;
	void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount) override;

	// NT traversal requires force exchange!
	bool requiresForceExchange() const override { return true; }

	// NT has no upper limit for number of cells in cutoff!
	unsigned maxCellsInCutoff() const override { return std::numeric_limits<unsigned>::max(); }

protected:
	void processBaseCell(CellProcessor& cellProcessor, unsigned long baseIndex) const;

	// All pairs that have to be processed when calculating the forces (excluding self)
	std::array<std::pair<long, long>, 62> _cellPairOffsets;
	std::array<std::pair<std::array<long, 3>, std::array<long, 3>>, 62> _offsets3D;

private:
	void computeOffsets();
	void computeOffsets3D(double cellLength[3], double cutoff);
};

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::computeOffsets3D(double cellLength[3], double cutoff) {}

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::computeOffsets() {}

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::traverseCellPairs(CellProcessor& cellProcessor) {
	unsigned long start = 0ul;
	unsigned long end = this->_cells->size();
	for (auto i = start; i < end; ++i) {
		processBaseCell(cellProcessor, i);
	}
}

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::traverseCellPairsOuter(CellProcessor& cellProcessor) {
	global_log->error() << "NT: overlapping Comm not implemented." << std::endl;
	Simulation::exit(46);
}

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage,
																	 unsigned stageCount) {
	global_log->error() << "NT: overlapping Comm not implemented." << std::endl;
	Simulation::exit(47);
}

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::processBaseCell(CellProcessor& cellProcessor,
															  unsigned long baseIndex) const {
	unsigned long maxIndex = this->_cells->size() - 1;

	CellTemplate& baseCell = this->_cells->at(baseIndex);

	if (not baseCell.isHaloCell()) {
		// Process all cell pairs for this cell
		for (auto& current_pair : _cellPairOffsets) {
			unsigned long offset1 = current_pair.first;
			unsigned long cellIndex1 = baseIndex + offset1;

			unsigned long offset2 = current_pair.second;
			unsigned long cellIndex2 = baseIndex + offset2;

			CellTemplate& cell1 = this->_cells->at(cellIndex1);
			CellTemplate& cell2 = this->_cells->at(cellIndex2);

			const bool sumAllMacroscopic = true;
			cellProcessor.processCellPair(cell1, cell2, sumAllMacroscopic);
		}

		// Process base cell itself
		cellProcessor.processCell(baseCell);
	}
}
