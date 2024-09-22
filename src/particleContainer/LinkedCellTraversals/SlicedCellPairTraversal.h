/*
 * SlicedCellPairTraversal.h
 *
 *  Created on: 16 May 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_SLICEDCELLPAIRTRAVERSAL_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_SLICEDCELLPAIRTRAVERSAL_H_

#include "particleContainer/LinkedCellTraversals/C08BasedTraversals.h"
#include "utils/ThreeElementPermutations.h"
#include "WrapOpenMP.h"

struct SlicedCellPairTraversalData : CellPairTraversalData {
};

template<class CellTemplate>
class SlicedCellPairTraversal: public C08BasedTraversals<CellTemplate> {
public:
	SlicedCellPairTraversal(
		std::vector<CellTemplate> &cells, const std::array<unsigned long, 3> &dims) :
		C08BasedTraversals<CellTemplate>(cells, dims)
		, _locks(mardyn_get_max_threads() - 1, nullptr)
	{
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			int myid = mardyn_get_thread_num();
			if(myid != mardyn_get_max_threads()-1) {
				_locks[myid] = new mardyn_lock_t();
			}
		}
	}

	using C08BasedTraversals<CellTemplate>::rebuild;

	virtual ~SlicedCellPairTraversal() {
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			int myid = mardyn_get_thread_num();
			if(myid != mardyn_get_max_threads()-1) {
				delete _locks[myid];
			}
		}
	}

	void traverseCellPairs(CellProcessor& cellProcessor);
	void traverseCellPairsOuter(CellProcessor& cellProcessor);
	void traverseCellPairsInner(CellProcessor& cellProcessor, unsigned stage, unsigned stageCount);

	static bool isApplicable(const std::array<unsigned long, 3>& dims);

	static bool isApplicable(
		const std::array<unsigned long, 3>& start,
		const std::array<unsigned long, 3>& end);

protected:
	void traverseCellPairsBackend(CellProcessor& cellProcessor,
		const std::array<unsigned long, 3>& start,
		const std::array<unsigned long, 3>& end
		);

private:
	enum LockType {
		MY_LOCK,
		NEXT_LOCK
	};

	void releaseMyLock();
	void acquireLock(LockType l);
	void initLocks();
	void destroyLocks();

	std::vector<mardyn_lock_t *> _locks;
};

template<class CellTemplate>
inline void SlicedCellPairTraversal<CellTemplate>::traverseCellPairs(
		CellProcessor& cellProcessor) {
	using std::array;
	const std::array<unsigned long, 3> start = { 0, 0, 0 };
	std::array<unsigned long, 3> end;
	for (int d = 0; d < 3; ++d) {
		end[d] = this->_dims[d] - 1;
	}

	traverseCellPairsBackend(cellProcessor, start, end);
}

template<class CellTemplate>
inline void SlicedCellPairTraversal<CellTemplate>::traverseCellPairsOuter(
		CellProcessor& cellProcessor) {
	using std::array;

	{
		unsigned long minsize = std::min(this->_dims[0], std::min(this->_dims[1], this->_dims[2]));

		if (minsize <= 5) {
			// iterating in the inner region didn't do anything. Iterate normally.
			traverseCellPairs(cellProcessor);
			return;
		}
	}

	// values, which are modified by 2 are actually modified by the respective stride, which is always 2
	// lower part
	std::array<unsigned long, 3> startZlo = { 0, 0, 0 };
	std::array<unsigned long, 3> endZlo = { this->_dims[0] - 1, this->_dims[1] - 1, 2 };
	traverseCellPairsBackend(cellProcessor, startZlo, endZlo);

	// upper part
	std::array<unsigned long, 3> startZhi = { 0, 0, this->_dims[2] - 3 };
	std::array<unsigned long, 3> endZhi = { this->_dims[0] - 1, this->_dims[1] - 1, this->_dims[2] - 1 };
	traverseCellPairsBackend(cellProcessor, startZhi, endZhi);

	// halo & boundaries in y direction
	// boundaries in z direction are excluded!
	// lower part
	std::array<unsigned long, 3> startYlo = { 0, 0, 2 };
	std::array<unsigned long, 3> endYlo = { this->_dims[0] - 1, 2, this->_dims[2] - 3 };
	traverseCellPairsBackend(cellProcessor, startYlo, endYlo);

	// upper part
	std::array<unsigned long, 3> startYhi = { 0, this->_dims[1] - 3, 2 };
	std::array<unsigned long, 3> endYhi = { this->_dims[0] - 1, this->_dims[1] - 1, this->_dims[2] - 3 };
	traverseCellPairsBackend(cellProcessor, startYhi, endYhi);

	// halo & boundaries in x direction
	// boundaries in z and y direction are excluded!
	// lower part
	std::array<unsigned long, 3> startXlo = { 0, 2, 2 };
	std::array<unsigned long, 3> endXlo = { 2, this->_dims[1] - 3, this->_dims[2] - 3 };
	traverseCellPairsBackend(cellProcessor, startXlo, endXlo);

	// upper part
	std::array<unsigned long, 3> startXhi = { this->_dims[0] - 3, 2, 2 };
	std::array<unsigned long, 3> endXhi = { this->_dims[0] - 1, this->_dims[1] - 3, this->_dims[2] - 3 };
	traverseCellPairsBackend(cellProcessor, startXhi, endXhi);
}

template<class CellTemplate>
inline void SlicedCellPairTraversal<CellTemplate>::traverseCellPairsInner(
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

	traverseCellPairsBackend(cellProcessor, lower, upper);
}

template<class CellTemplate>
inline void SlicedCellPairTraversal<CellTemplate>::traverseCellPairsBackend(
		CellProcessor& cellProcessor,
		const std::array<unsigned long, 3>& start,
		const std::array<unsigned long, 3>& end
		) {

	using namespace Permute3Elements;
	using std::array;

	// Note: in the following we quasi-reimplement an OpenMP for-loop parallelisation with static scheduling
	if (not isApplicable(start, end) ) {
		Log::global_log->error() << "The SlicedCellPairTraversal is not applicable. Aborting." << std::endl;
		MARDYN_EXIT(1);
	}

	std::array<unsigned long, 3> diff;

	unsigned long num_cells = 1;
	for(int d = 0; d < 3; ++d) {
		mardyn_assert(end[d] >= start[d]);
		unsigned long difference = end[d] - start[d];
		num_cells *= difference;
		diff[d] = difference;
	}

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		Permutation perm = getPermutationForIncreasingSorting(diff);
		std::array<unsigned long, 3> diff_permuted = permuteForward(perm, diff);
		std::array<unsigned long, 3> start_permuted = permuteForward(perm, start);

		initLocks();

		int my_id = mardyn_get_thread_num();
		int num_threads = mardyn_get_num_threads();

		const unsigned long my_start = num_cells * my_id / num_threads;
		const unsigned long my_end = num_cells * (my_id + 1) / num_threads;
		const unsigned long my_num_cells = my_end - my_start; // a rough measure should be enough?
		const unsigned long slice_size = diff_permuted[0] * diff_permuted[1];
		unsigned long my_progress_counter = 0;

		acquireLock(MY_LOCK);

		#if defined(_OPENMP)
		#pragma omp barrier
		#endif

		// Note: manually implementing omp for schedule(static) collapse(3) unfortunately
		// because I am relying on the execution order, i.e. thread 0 processes first chunk, thread 1 processes second chunk
		// and so on, which I don't know whether is guaranteed by the OpenMP standard.
		for (unsigned long i = my_start; i < my_end; ++i) {
			if (my_progress_counter == my_num_cells - slice_size) {
				acquireLock(NEXT_LOCK);
			}

			// unroll
			std::array<unsigned long, 3> newInd_permuted = threeDimensionalMapping::oneToThreeD(i,diff_permuted);

			// add inner-, middle- and outer-start
			for(int d = 0; d < 3; ++d) {
				newInd_permuted[d] += start_permuted[d];
			}

			// permute newInd backwards
			std::array<unsigned long, 3> newInd = permuteBackward(perm, newInd_permuted);

			// get actual index
			unsigned long cellIndex = threeDimensionalMapping::threeToOneD(newInd, this->_dims);

			this->processBaseCell(cellProcessor, cellIndex);

			++my_progress_counter;

			if (my_progress_counter == slice_size) {
				releaseMyLock();
			}
		}

		#if defined(_OPENMP)
		#pragma omp barrier
		#endif

		destroyLocks();
	} /* omp parallel */
}

template<class CellTemplate>
inline void SlicedCellPairTraversal<CellTemplate>::releaseMyLock() {
//	#if defined(_OPENMP)
//	#pragma omp parallel
//	#endif
	{
		int myid = mardyn_get_thread_num();
		if(myid != 0) {
			mardyn_unset_lock(_locks[myid - 1]);
		}
	}
}

template<class CellTemplate>
inline void SlicedCellPairTraversal<CellTemplate>::acquireLock(LockType l) {
//	#if defined(_OPENMP)
//	#pragma omp parallel
//	#endif
	{
		int myid = mardyn_get_thread_num();
		switch(l) {
		case MY_LOCK:
			if (myid != 0) {
				mardyn_set_lock(_locks[myid - 1]);
			}
			break;
		case NEXT_LOCK:
			if (myid != mardyn_get_num_threads()-1) {
				mardyn_set_lock(_locks[myid]);
			}
			break;
		}
	}
}

template<class CellTemplate>
inline void SlicedCellPairTraversal<CellTemplate>::initLocks() {
//	#if defined(_OPENMP)
//	#pragma omp parallel
//	#endif
	{
		int myid = mardyn_get_thread_num();
		if(myid != 0) {
			mardyn_init_lock(_locks[myid - 1]);
		}
	}
}

template<class CellTemplate>
inline bool SlicedCellPairTraversal<CellTemplate>::isApplicable(
		const std::array<unsigned long, 3>& dims) {
	using std::array;
	const std::array<unsigned long, 3> start = { 0, 0, 0 };
	std::array<unsigned long, 3> end;
	for (int d = 0; d < 3; ++d) {
		end[d] = dims[d] - 1;
	}

	return isApplicable(start, end);
}

template<class CellTemplate>
inline bool SlicedCellPairTraversal<CellTemplate>::isApplicable(
		const std::array<unsigned long, 3>& start,
		const std::array<unsigned long, 3>& end) {

	using namespace Permute3Elements;
	using std::array;

	std::array<unsigned long, 3> dims = {end[0] - start[0], end[1] - start[1], end[2] - start[2]};
	Permutation permutation = getPermutationForIncreasingSorting(dims);
	std::array<unsigned long, 3> dimsPermuted = permuteForward(permutation, dims);

	int num_threads = mardyn_get_max_threads();

	size_t num_cells = 1;
	for(int d = 0; d < 3; ++d) {
		num_cells *= end[d] - start[d];
	}

	const size_t my_num_cells = num_cells / num_threads; // a rough measure should be enough? // lower bound ?
	const size_t slice_size = dimsPermuted[0] * dimsPermuted[1];

	const bool ret = my_num_cells >= 2 * slice_size or num_threads == 1;

	return ret;
}

template<class CellTemplate>
inline void SlicedCellPairTraversal<CellTemplate>::destroyLocks() {
//	#if defined(_OPENMP)
//	#pragma omp parallel
//	#endif
	{
		int myid = mardyn_get_thread_num();
		if(myid != mardyn_get_num_threads()-1) {
			mardyn_destroy_lock(_locks[myid]);
		}
	}
}

#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_SLICEDCELLPAIRTRAVERSAL_H_ */
