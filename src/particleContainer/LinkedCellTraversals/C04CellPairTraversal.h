/*
 * C04CellPairTraversal.h
 *
 *  Created on: 26 Mar 2018
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C04CELLPAIRTRAVERSAL_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C04CELLPAIRTRAVERSAL_H_

#include "C08BasedTraversals.h"
#include "utils/threeDimensionalMapping.h"
#include "utils/mardyn_assert.h"

struct C04CellPairTraversalData : CellPairTraversalData {
};

template <class CellTemplate>
class C04CellPairTraversal: public C08BasedTraversals<CellTemplate> {
public:
	C04CellPairTraversal(std::vector<CellTemplate>& cells,
			const std::array<unsigned long, 3>& dims) :
			C08BasedTraversals<CellTemplate>(cells, dims) {
		computeOffsets32Pack();
	}
	virtual ~C04CellPairTraversal() {
	}

	virtual void rebuild(std::vector<CellTemplate> &cells,
			 const std::array<unsigned long, 3> &dims, double cellLength[3], double cutoff,
			 CellPairTraversalData *data) {
		C08BasedTraversals<CellTemplate>::rebuild(cells, dims, cellLength, cutoff, data);
		computeOffsets32Pack();
	}

	void traverseCellPairs(CellProcessor& cellProcessor);
	void traverseCellPairsOuter(CellProcessor& cellProcessor) {}
	void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount) {}

private:
	void traverseCellPairsBackend(CellProcessor& cellProcessor,
			const std::array<long, 3>& start,
			const std::array<long, 3>& end);

	void traverseSingleColor(CellProcessor& cellProcessor,
			int color,
			const std::array<long, 3>& start,
			const std::array<long, 3>& end);

	void processBasePack32(CellProcessor& cellProcessor,
			const std::array<long, 3>& base3DIndex,
			const std::array<long, 3>& start,
			const std::array<long, 3>& end) const;

	void computeOffsets32Pack();

	long parity(long x, long y, long z) const {
		return (x + y + z + 24) % 8;
	}


	std::array<std::array<long, 3>, 32> _cellOffsets32Pack;
};

template<class CellTemplate>
inline void C04CellPairTraversal<CellTemplate>::traverseCellPairs(
		CellProcessor& cellProcessor) {
	std::array<long, 3> start, end;
	for (int d = 0; d < 3; ++d) {
		start[d] = 0l;
		end[d] = static_cast<long>(this->_dims[d]) - 1;
	}
	traverseCellPairsBackend(cellProcessor, start, end);
}

template<class CellTemplate>
void C04CellPairTraversal<CellTemplate>::computeOffsets32Pack() {
	using threeDimensionalMapping::threeToOneD;
	using std::make_pair;

	int i = 0;
	long z = 0l;
	_cellOffsets32Pack[i++] = {1l, 1l, z};
	_cellOffsets32Pack[i++] = {1l, 2l, z};
	_cellOffsets32Pack[i++] = {2l, 1l, z};
	_cellOffsets32Pack[i++] = {2l, 2l, z};

	/* z = 1ul; z = 2ul */
	for (z = 1l; z < 3l; ++z) {
		for (long y = 0l; y < 4l; y++) {
			for (long x = 0l; x < 4l; x++) {
				if ((x == 0l and y == 0l) or
					(x == 3l and y == 0l) or
					(x == 0l and y == 3l) or
					(x == 3l and y == 3l)) {
					continue;
				}
				_cellOffsets32Pack[i++] = {x, y, z};
			}
		}
	}


	z = 3ul;
	_cellOffsets32Pack[i++] = {1l, 1l, z};
	_cellOffsets32Pack[i++] = {1l, 2l, z};
	_cellOffsets32Pack[i++] = {2l, 1l, z};
	_cellOffsets32Pack[i++] = {2l, 2l, z};

	mardyn_assert(i == 32);

}

template<class CellTemplate>
void C04CellPairTraversal<CellTemplate>::processBasePack32(
		CellProcessor& cellProcessor,
		const std::array<long, 3>& base3DIndex,
		const std::array<long, 3>& start,
		const std::array<long, 3>& end
		) const {
	using threeDimensionalMapping::threeToOneD;
	std::array<long, 3> index;
	std::array<long, 3> signedDims;
	for (int d = 0; d < 3; ++d) {
		signedDims[d] = static_cast<long>(this->_dims[d]);
	}

	for (auto Offset32Pack : _cellOffsets32Pack) {
		// compute 3D index
		bool isIn = true;
		for (int d = 0; d < 3; ++d) {
			index[d] = base3DIndex[d] + Offset32Pack[d];
			isIn = isIn and (index[d] >= start[d]) and (index[d] < end[d]);
		}

		if (isIn) {
			unsigned long ulIndex = static_cast<unsigned long>(threeToOneD(index, signedDims));
			C08BasedTraversals<CellTemplate>::processBaseCell(cellProcessor, ulIndex);
		} else {
			continue;
		}
	}
}

template<class CellTemplate>
void C04CellPairTraversal<CellTemplate>::traverseCellPairsBackend(
		CellProcessor & cellProcessor,
		const std::array<long, 3>& start,
		const std::array<long, 3>& end) {

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{

		for (int color = 0; color < 4; ++color) {

			traverseSingleColor(cellProcessor, color, start, end);

			#if defined(_OPENMP)
			if (color < 3) {
				#pragma omp barrier
			}
			#endif
		}
	} /* close parallel region */
}

template<class CellTemplate>
void C04CellPairTraversal<CellTemplate>::traverseSingleColor(CellProcessor& cellProcessor,
		int color,
		const std::array<long, 3>& start,
		const std::array<long, 3>& end) {

	std::array<long, 3> intersectionStart;
	for (int d = 0; d < 3; ++d) {
		intersectionStart[d] = start[d] - 2;
	}

	// we need to traverse one BCC grid, which consists of two cartesian grids

	// colors 0 and 2 form one cartesian grid
	// colors 1 and 3 form another cartesian grid, whose origin is shifted by (2,2,2)

	// determine a starting point of one of the grids
	std::array<long, 3> startOfThisColor {0l, 0l, 0l};

	switch(color % 2) {
	case 0:
		// colours 0 and 2
		startOfThisColor = intersectionStart;
		break;
	case 1:
		// colours 1 and 3
		startOfThisColor = start;
		break;
	default:
		mardyn_assert(false);
	}

	long correctParity;
	correctParity = parity(startOfThisColor[0], startOfThisColor[1], startOfThisColor[2]);
	if (color >= 2) {
		correctParity += 4;
	}

	// to fix compiler complaints about perfectly nested loop.
	const long startX = startOfThisColor[0], endX = end[0];
	const long startY = startOfThisColor[1], endY = end[1];
	const long startZ = startOfThisColor[2], endZ = end[2];

	const auto loop_size = static_cast<size_t>(std::ceil(static_cast<double>(endX - startX) / 4) *
											   std::ceil(static_cast<double>(endY - startY) / 4) *
											   std::ceil(static_cast<double>(endZ - startZ) / 4));

	// Here, we use a smaller max_chunk_size compared to C08, as c04 work items are bigger.
	const int chunk_size = chunk_size::getChunkSize(loop_size, 10000, 20);

	// first cartesian grid
	#if defined(_OPENMP)
	#pragma omp for schedule(dynamic, chunk_size) collapse(3) nowait
	#endif
	for (long z = startZ; z < endZ; z += 4) {
		for (long y = startY; y < endY; y += 4) {
			for (long x = startX; x < endX; x += 4) {

				long par = parity(x, y, z);

				if (par != correctParity) {
					continue;
				}

				std::array<long, 3> base3DIndex = {x, y, z};
				processBasePack32(cellProcessor, base3DIndex, start, end);
			}
		}
	}
}

#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C04CELLPAIRTRAVERSAL_H_ */
