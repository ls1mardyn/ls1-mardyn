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

/**
 * This is the traversal corresponding to the neutral territory zonal method scheme.
 * Interactions occur between cells of the plate and cells of the tower.
 * Algorithm description in:
 * D. E. Shaw: A Fast, Scalable Method for the Parallel Evaluation of Distance-Limited Pairwise Particle Interactions.
 * Description of implementation in:
 * S. Seckler: Algorithm and Performance Engineering for HPC Particle Simulations
 * @tparam CellTemplate
 */
template <class CellTemplate>
class NeutralTerritoryTraversal : public CellPairTraversals<CellTemplate> {
public:
	NeutralTerritoryTraversal(std::vector<CellTemplate>& cells, const std::array<unsigned long, 3>& dims,
							  double cellLength[3], double cutoff)
		: CellPairTraversals<CellTemplate>(cells, dims) {
		computeOffsets3D(cellLength, cutoff);
		computeOffsets();
	}
	~NeutralTerritoryTraversal() override = default;

	/**
	 * Reset all necessary data.
	 */
	void rebuild(std::vector<CellTemplate>& cells, const std::array<unsigned long, 3>& dims, double cellLength[3],
				 double cutoff, CellPairTraversalData* data) override {
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
	std::vector<std::pair<long, long>> _cellPairOffsets;
	std::vector<std::pair<std::array<long, 3>, std::array<long, 3>>> _offsets3D;

private:
	void computeOffsets();
	void computeOffsets3D(const double cellLength[3], double cutoff);
};

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::computeOffsets3D(const double cellLength[3], const double cutoff) {
	_offsets3D.clear();
	_cellPairOffsets.clear();
	long num_x = std::ceil(cutoff / cellLength[0]);
	long num_y = std::ceil(cutoff / cellLength[1]);
	long num_z = std::ceil(cutoff / cellLength[2]);
	for (long plate_x = 0; plate_x <= num_x; ++plate_x) {
		// Start plate_y at 0, iff plate_x == 0.
		long start_plate_y = plate_x == 0 ? 0 : -num_y;
		for (long plate_y = start_plate_y; plate_y <= num_y; ++plate_y) {
			// Start tower_z at 1, iff plate_x == 0 and plate_y == 0.
			// This also prevents both tower and plate to be {0,0,0}.
			long start_tower_z = (plate_x == 0 and plate_y == 0) ? 1 : -num_z;
			for (long tower_z = start_tower_z; tower_z <= num_z; ++tower_z) {
				std::array<long, 3> plate{plate_x, plate_y, 0};
				std::array<long, 3> tower{0, 0, tower_z};
				_offsets3D.emplace_back(tower, plate);
			}
		}
	}
}

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::computeOffsets() {
	// Clear _cellPairOffsets!
	_cellPairOffsets.clear();

	using threeDimensionalMapping::threeToOneD;

	// Dim array is int but we need it as long for some reason (copied from C08BasedTraversal)
	std::array<long, 3> dims{};
	for (int d = 0; d < 3; ++d) {
		dims[d] = static_cast<long>(this->_dims[d]);
	}

	for (const auto& [firstOffset, secondOffset] : _offsets3D) {
		auto [firstOffsetX, firstOffsetY, firstOffsetZ] = firstOffset;
		auto [secondOffsetX, secondOffsetY, secondOffsetZ] = secondOffset;

		// convert 3d index to 1d
		auto aIndex = threeToOneD(firstOffsetX, firstOffsetY, firstOffsetZ, dims);
		auto bIndex = threeToOneD(secondOffsetX, secondOffsetY, secondOffsetZ, dims);

		// store offset pair
		_cellPairOffsets.emplace_back(aIndex, bIndex);
	}
}

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
	std::ostringstream error_message;
	error_message << "NT: overlapping Comm not implemented." << std::endl;
	MARDYN_EXIT(error_message.str());
}

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage,
																	 unsigned stageCount) {
	std::ostringstream error_message;
	error_message << "NT: overlapping Comm not implemented." << std::endl;
	MARDYN_EXIT(error_message.str());
}

template <class CellTemplate>
void NeutralTerritoryTraversal<CellTemplate>::processBaseCell(CellProcessor& cellProcessor,
															  unsigned long baseIndex) const {
	CellTemplate& baseCell = this->_cells->at(baseIndex);

	if (not baseCell.isHaloCell()) {
		// Process all cell pairs for this cell
		for (const auto& [offset1, offset2] : _cellPairOffsets) {
			unsigned long cellIndex1 = baseIndex + offset1;
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
