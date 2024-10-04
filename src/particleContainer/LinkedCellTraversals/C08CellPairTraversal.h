/*
 * C08CellPairTraversal.h
 *
 *  Created on: 15 May 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08CELLPAIRTRAVERSAL_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08CELLPAIRTRAVERSAL_H_

#include "C08BasedTraversals.h"
#include "utils/GetChunkSize.h"
#include "utils/mardyn_assert.h"
#include "utils/threeDimensionalMapping.h"

struct C08CellPairTraversalData : CellPairTraversalData {
};

template <class CellTemplate, bool eighthShell = false>
class C08CellPairTraversal : public C08BasedTraversals<CellTemplate> {
public:
	C08CellPairTraversal(
			std::vector<CellTemplate>& cells,
			const std::array<unsigned long, 3>& dims) :
			C08BasedTraversals<CellTemplate>(cells, dims) {
	}
	~C08CellPairTraversal() = default;

	using C08BasedTraversals<CellTemplate>::rebuild;

	void traverseCellPairs(CellProcessor& cellProcessor) override;
	void traverseCellPairsOuter(CellProcessor& cellProcessor) override;
	void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount) override;

	bool requiresForceExchange() const override {return eighthShell;}

private:
	void traverseCellPairsBackend(CellProcessor& cellProcessor,
			const std::array<unsigned long, 3> & start,
			const std::array<unsigned long, 3> & end,
			const std::array<unsigned long, 3> & stride) const;
};


template<class CellTemplate, bool eighthShell>
void C08CellPairTraversal<CellTemplate, eighthShell>::traverseCellPairs(
		CellProcessor& cellProcessor) {

	using std::array;
	const std::array<unsigned long, 3> strides = { 2, 2, 2 };
	std::array<unsigned long, 3> end;
	for (int d = 0; d < 3; ++d) {
		end[d] = this->_dims[d] - 1;
	}

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		for (unsigned long col = 0; col < 8; ++col) {
			std::array<unsigned long, 3> begin = threeDimensionalMapping::oneToThreeD(col, strides);
			if (eighthShell) {
				// if we are using eighth shell, we start at 1,1,1 instead of 0,0,0
				for (unsigned short i = 0; i < 3; ++i) {
					begin[i] += 1;
				}
			}
			traverseCellPairsBackend(cellProcessor, begin, end, strides);
			#if defined(_OPENMP)
			#pragma omp barrier
			#endif
			// this barrier is needed, since we have a nowait in the backend
		}
	}
}

template<class CellTemplate, bool eighthShell>
void C08CellPairTraversal<CellTemplate, eighthShell>::traverseCellPairsOuter(
		CellProcessor& cellProcessor) {
	if(eighthShell){
		std::ostringstream error_message;
		error_message << "eightshell + overlapping not yet supported." << std::endl;
		MARDYN_EXIT(error_message);
	}
	using std::array;

	{
		unsigned long minsize = std::min(this->_dims[0], std::min(this->_dims[1], this->_dims[2]));

		if (minsize <= 5) {
			// iterating in the inner region didn't do anything. Iterate normally.
			traverseCellPairs(cellProcessor);
			return;
		}
	}

	const std::array<unsigned long, 3> strides2 = { 2, 2, 2 };

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		for (unsigned long col = 0; col < 8; ++col) {
			// halo & boundaries in z direction
			const std::array< unsigned long, 3> begin = threeDimensionalMapping::oneToThreeD(col, strides2);

			// values, which are modified by 2 are actually modified by the respective stride, which is always 2
			std::array<unsigned long, 3> startZ = { begin[0], begin[1], begin[2] }; 									// default
			std::array<unsigned long, 3> endZ = { this->_dims[0] - 1, this->_dims[1] - 1, this->_dims[2] - 1 };		// default
			std::array<unsigned long, 3> stridesZ = {strides2[0], strides2[1], this->_dims[2] - 3};					// mod z
			traverseCellPairsBackend(cellProcessor, startZ, endZ, stridesZ);

			// halo & boundaries in y direction
			// boundaries in z direction are excluded!
			std::array<unsigned long, 3> startY = { begin[0], begin[1], begin[2] + 2 };								// mod z
			std::array<unsigned long, 3> endY = { this->_dims[0] - 1, this->_dims[1] - 1, this->_dims[2] - 3 };		// mod z
			std::array<unsigned long, 3> stridesY = {strides2[0], this->_dims[1] - 3, strides2[2]};					// mod y
			traverseCellPairsBackend(cellProcessor, startY, endY, stridesY);

			// halo & boundaries in x direction
			// boundaries in z and y direction are excluded!
			std::array<unsigned long, 3> startX = { begin[0], begin[1] + 2, begin[2] + 2 };							// mod yz
			std::array<unsigned long, 3> endX = { this->_dims[0] - 1, this->_dims[1] - 3, this->_dims[2] - 3 };		// mod yz
			std::array<unsigned long, 3> stridesX = {this->_dims[0] - 3, strides2[1], strides2[2]};					// mod x
			traverseCellPairsBackend(cellProcessor, startX, endX, stridesX);

			#if defined(_OPENMP)
				// this barrier is needed, since we have a nowait in the backend
				// except at the last run - there we have the barrier at the end of the parallel region
				if (col < 7) {
					#pragma omp barrier
				}
			#endif
		}
	} // end pragma omp parallel
}

template<class CellTemplate, bool eighthShell>
void C08CellPairTraversal<CellTemplate, eighthShell>::traverseCellPairsInner(
		CellProcessor& cellProcessor, unsigned stage,
		unsigned stageCount) {
	using std::array;

	unsigned long splitdim = 0;
	unsigned long maxcellsize = this->_dims[0];
	for (unsigned long i = 1; i < 3; i++){
		if (this->_dims[i] > maxcellsize) {
			splitdim = i;
			maxcellsize = this->_dims[i];
		}
	}
	unsigned long splitsize = maxcellsize - 5;
	unsigned long minsize = std::min(this->_dims[0], std::min(this->_dims[1], this->_dims[2]));

	mardyn_assert(minsize >= 4);  // there should be at least 4 cells in each dimension, otherwise we did something stupid!

	if (minsize <= 5) {
		return;  // we can not iterate over any inner cells, that do not depend on boundary or halo cells
	}

	std::array<unsigned long, 3> lower;
	std::array<unsigned long, 3> upper;
	for (unsigned long i = 0; i < 3; i++) {
		lower[i] = 2;
		upper[i] = this->_dims[i] - 3;
	}
	lower[splitdim] = 2 + splitsize * stage / stageCount;  // at least 2
	upper[splitdim] = 2 + splitsize * (stage+1) / stageCount;  // at most _cellsPerDimension[i] - 3


	#if defined(_OPENMP)
		#pragma omp parallel
	#endif
	{
		std::array<unsigned long, 3> strides = {2, 2, 2};

		for (unsigned long col = 0; col < 8; ++col) {
			std::array<unsigned long, 3> startIndices = threeDimensionalMapping::oneToThreeD(col, strides);
			for (int i = 0; i < 3; i++) {
				startIndices[i] = startIndices[i] + lower[i];
			}

			traverseCellPairsBackend(cellProcessor, startIndices, upper, strides);
			#if defined(_OPENMP)
				// this barrier is needed, since we have a nowait in the backend
				// except at the last run - there we have the barrier at the end of the parallel region
				if (col < 7) {
					#pragma omp barrier
				}
			#endif
		}
	} // end pragma omp parallel
}


template<class CellTemplate, bool eighthShell>
void C08CellPairTraversal<CellTemplate, eighthShell>::traverseCellPairsBackend(CellProcessor& cellProcessor,
		const std::array<unsigned long, 3>& start, const std::array<unsigned long, 3>& end,
		const std::array<unsigned long, 3>& stride) const {

	// note parallel region is open outside

	// intel compiler demands following:
	// note parallel region is open outside
	const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
	const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];
	const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

	// number of iterations:
	const auto loop_size = static_cast<size_t>(std::ceil(static_cast<double>(end_x - start_x) / stride_x) *
											   std::ceil(static_cast<double>(end_y - start_y) / stride_y) *
											   std::ceil(static_cast<double>(end_z - start_z) / stride_z));

	// magic numbers: empirically determined to be somewhat efficient.
	const int chunk_size = chunk_size::getChunkSize(loop_size, 10000, 100);

#if defined(_OPENMP)
	#pragma omp for schedule(dynamic, chunk_size) collapse(3) nowait
	#endif
	for (unsigned long z = start_z; z < end_z; z += stride_z) {
		for (unsigned long y = start_y; y < end_y; y += stride_y) {
			for (unsigned long x = start_x; x < end_x; x += stride_x) {
				unsigned long baseIndex = threeDimensionalMapping::threeToOneD(x, y, z, this->_dims);
				C08BasedTraversals<CellTemplate>::template processBaseCell<eighthShell>(cellProcessor, baseIndex);
			}
		}
	}
}


#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_C08CELLPAIRTRAVERSAL_H_ */
