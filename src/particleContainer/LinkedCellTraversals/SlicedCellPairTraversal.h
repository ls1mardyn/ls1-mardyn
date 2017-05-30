/*
 * SlicedCellPairTraversal.h
 *
 *  Created on: 16 May 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_SLICEDCELLPAIRTRAVERSAL_H_
#define SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_SLICEDCELLPAIRTRAVERSAL_H_

#include "particleContainer/LinkedCellTraversals/C08BasedTraversals.h"
#include "WrapOpenMP.h"

template<class CellTemplate>
class SlicedCellPairTraversal: public C08BasedTraversals<CellTemplate> {
public:
	SlicedCellPairTraversal(
		std::vector<CellTemplate>& cells, std::array<unsigned long, 3>& dims) :
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

	enum Permutation {
		XYZ,
		XZY,
		YZX,
		YXZ,
		ZXY,
		ZYX
	};
	static void arrangeDimensions(
		unsigned long & outer_start, unsigned long & middle_start, unsigned long & inner_start,
		unsigned long & outer_end, unsigned long & middle_end, unsigned long & inner_end,
		Permutation & permutation,
		const std::array<unsigned long, 3>& start,
		const std::array<unsigned long, 3>& end);

	unsigned long mapCellIndex(
		unsigned long inner, unsigned long middle, unsigned long outer,
		Permutation perm) const;

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
	const array<unsigned long, 3> start = { 0, 0, 0 };
	array<unsigned long, 3> end;
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
		unsigned long minsize = min(this->_dims[0], min(this->_dims[1], this->_dims[2]));

		if (minsize <= 5) {
			// iterating in the inner region didn't do anything. Iterate normally.
			traverseCellPairs(cellProcessor);
			return;
		}
	}

	// values, which are modified by 2 are actually modified by the respective stride, which is always 2
	// lower part
	array<unsigned long, 3> startZlo = { 0, 0, 0 };
	array<unsigned long, 3> endZlo = { this->_dims[0] - 1, this->_dims[1] - 1, 2 };
	traverseCellPairsBackend(cellProcessor, startZlo, endZlo);

	// upper part
	array<unsigned long, 3> startZhi = { 0, 0, this->_dims[2] - 3 };
	array<unsigned long, 3> endZhi = { this->_dims[0] - 1, this->_dims[1] - 1, this->_dims[2] - 1 };
	traverseCellPairsBackend(cellProcessor, startZhi, endZhi);

	// halo & boundaries in y direction
	// boundaries in z direction are excluded!
	// lower part
	array<unsigned long, 3> startYlo = { 0, 0, 2 };
	array<unsigned long, 3> endYlo = { this->_dims[0] - 1, 2, this->_dims[2] - 3 };
	traverseCellPairsBackend(cellProcessor, startYlo, endYlo);

	// upper part
	array<unsigned long, 3> startYhi = { 0, this->_dims[1] - 3, 2 };
	array<unsigned long, 3> endYhi = { this->_dims[0] - 1, this->_dims[1] - 1, this->_dims[2] - 3 };
	traverseCellPairsBackend(cellProcessor, startYhi, endYhi);

	// halo & boundaries in x direction
	// boundaries in z and y direction are excluded!
	// lower part
	array<unsigned long, 3> startXlo = { 0, 2, 2 };
	array<unsigned long, 3> endXlo = { 2, this->_dims[1] - 3, this->_dims[2] - 3 };
	traverseCellPairsBackend(cellProcessor, startXlo, endXlo);

	// upper part
	array<unsigned long, 3> startXhi = { this->_dims[0] - 3, 2, 2 };
	array<unsigned long, 3> endXhi = { this->_dims[0] - 1, this->_dims[1] - 3, this->_dims[2] - 3 };
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
	unsigned long minsize = min(this->_dims[0], min(this->_dims[1], this->_dims[2]));

	mardyn_assert(minsize >= 4);  // there should be at least 4 cells in each dimension, otherwise we did something stupid!

	if (minsize <= 5) {
		return;  // we can not iterate over any inner cells, that do not depend on boundary or halo cells
	}

	array<unsigned long, 3> lower;
	array<unsigned long, 3> upper;
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
	// Note: in the following we quasi-reimplement an OpenMP for-loop parallelisation with static scheduling
	mardyn_assert(isApplicable(start, end));

	size_t num_cells = 1;
	for(int d = 0; d < 3; ++d) {
		num_cells *= end[d] - start[d];
	}

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{

		unsigned long outer_start, middle_start, inner_start;
		unsigned long outer_end, middle_end, inner_end;
		Permutation permutation;

		// order dimensions from longest to shortest:
		arrangeDimensions(
			outer_start, middle_start, inner_start,
			outer_end, middle_end, inner_end,
			permutation,
			start, end);

		initLocks();

		int my_id = mardyn_get_thread_num();
		int num_threads = mardyn_get_num_threads();

		const size_t my_start = num_cells * my_id / num_threads;
		const size_t my_end = num_cells * (my_id + 1) / num_threads;
		const size_t my_num_cells = my_end - my_start; // a rough measure should be enough?
		const size_t slice_size = (middle_end - middle_start) * (inner_end - inner_start);
		size_t my_progress_counter = 0;

		acquireLock(MY_LOCK);

		#if defined(_OPENMP)
		#pragma omp barrier
		#pragma omp for schedule(static) collapse(3)
		#endif
		for (unsigned long outer = outer_start; outer < outer_end; ++outer) {
			for (unsigned long middle = middle_start; middle < middle_end; ++middle) {
				for (unsigned long inner = inner_start; inner < inner_end; ++inner) {
					if (my_progress_counter == my_num_cells - slice_size) {
						acquireLock(NEXT_LOCK);
					}
					unsigned long cellIndex = mapCellIndex(inner, middle, outer, permutation);
					this->processBaseCell(cellProcessor, cellIndex);

					++my_progress_counter;

					if (my_progress_counter == slice_size) {
						releaseMyLock();
					}
				}
			}
		}

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
		if(myid != mardyn_get_num_threads()-1) {
			mardyn_unset_lock(_locks[myid]);
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
	const array<unsigned long, 3> start = { 0, 0, 0 };
	array<unsigned long, 3> end;
	for (int d = 0; d < 3; ++d) {
		end[d] = dims[d] - 1;
	}

	return isApplicable(start, end);
}

template<class CellTemplate>
inline bool SlicedCellPairTraversal<CellTemplate>::isApplicable(
		const std::array<unsigned long, 3>& start,
		const std::array<unsigned long, 3>& end) {
	unsigned long outer_start, middle_start, inner_start;
	unsigned long outer_end, middle_end, inner_end;
	Permutation permutation;

	// order dimensions from longest to shortest:
	arrangeDimensions(
		outer_start, middle_start, inner_start,
		outer_end, middle_end, inner_end,
		permutation,
		start, end);

	int num_threads = mardyn_get_max_threads();

	size_t num_cells = 1;
	for(int d = 0; d < 3; ++d) {
		num_cells *= end[d] - start[d];
	}

	const size_t my_num_cells = num_cells / num_threads; // a rough measure should be enough? // lower bound ?
	const size_t slice_size = (middle_end - middle_start) * (inner_end - inner_start);

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

template<class CellTemplate>
void SlicedCellPairTraversal<CellTemplate>::arrangeDimensions(
		unsigned long & outer_start, unsigned long & middle_start, unsigned long & inner_start,
		unsigned long & outer_end, unsigned long & middle_end, unsigned long & inner_end,
		Permutation & permutation,
		const std::array<unsigned long, 3>& start,
		const std::array<unsigned long, 3>& end) {

	// shortcuts:
	int z_dim = end[2] - start[2];
	int y_dim = end[1] - start[1];
	int x_dim = end[0] - start[0];
	int z_start = start[2];
	int y_start = start[1];
	int x_start = start[0];

	enum {
		X = 0,
		Y = 1,
		Z = 2,
		bad = 9999
	} outer = bad, middle = bad, inner = bad;


	if (y_dim < x_dim) {
		if (z_dim < x_dim) {
			if (z_dim < y_dim) {
				permutation = ZYX;
			} else {
				permutation = YZX;
			}
		} else {
			permutation = YXZ;
		}
	} else {
		if (z_dim < y_dim) {
			if (z_dim < x_dim) {
				permutation = ZXY;
			} else {
				permutation = XZY;
			}
		} else {
			permutation = XYZ;
		}
	}

	switch(permutation) {
	case XYZ:
		inner = X; middle = Y; outer = Z;
		break;
	case XZY:
		inner = X; middle = Z; outer = Y;
		break;
	case YXZ:
		inner = Y; middle = X; outer = Z;
		break;
	case YZX:
		inner = Y; middle = Z; outer = X;
		break;
	case ZXY:
		inner = Z; middle = X; outer = Y;
		break;
	case ZYX:
		inner = Z; middle = Y; outer = X;
		break;
	}

	outer_end = end[outer];
	middle_end = end[middle];
	inner_end = end[inner];

	outer_start = start[outer];
	middle_start = start[middle];
	inner_start = start[inner];


	mardyn_assert(outer_end - outer_start >= middle_end - middle_start);
	mardyn_assert(middle_end - middle_start >= inner_end - inner_start);
}

template<class CellTemplate>
unsigned long SlicedCellPairTraversal<CellTemplate>::mapCellIndex(
		unsigned long inner, unsigned long middle, unsigned long outer,
		Permutation perm) const {
	unsigned long X, Y, Z;

	switch(perm) {
	case XYZ:
		X = inner; Y = middle; Z = outer;
		break;
	case XZY:
		X = inner; Z = middle; Y = outer;
		break;
	case YXZ:
		Y = inner; X = middle; Z = outer;
		break;
	case YZX:
		Y = inner; Z = middle; X = outer;
		break;
	case ZXY:
		Z = inner; X = middle; Y = outer;
		break;
	case ZYX:
		Z = inner; Y = middle; X = outer;
		break;
	}

	unsigned long ret = threeDimensionalMapping::threeToOneD(X, Y, Z, this->_dims);

	return ret;
}


#endif /* SRC_PARTICLECONTAINER_LINKEDCELLTRAVERSALS_SLICEDCELLPAIRTRAVERSAL_H_ */
