/*
 * WrapOpenMP.h
 *
 *  Created on: 15 Sep 2016
 *      Author: tchipevn
 */

#ifndef SRC_WRAPOPENMP_H_
#define SRC_WRAPOPENMP_H_

/**
 * Provide non-OpenMP versions of the most common OpenMP function calls,
 * so that they don't have to be wrapped in #ifdef-s every time.
 *
 * Proper wrapper and renaming necessary, because of -fopenmp-simd handling of gcc.
 *
 * Extend when necessary.
 */

#if defined(_OPENMP)
	#include <omp.h>

	inline int mardyn_get_thread_num()  {return omp_get_thread_num(); }
	inline int mardyn_get_num_threads() {return omp_get_num_threads();}
	inline int mardyn_get_max_threads() {return omp_get_max_threads();}
#else
	inline int mardyn_get_thread_num()  {return 0;}
	inline int mardyn_get_num_threads() {return 1;}
	inline int mardyn_get_max_threads() {return 1;}
#endif

#endif /* SRC_WRAPOPENMP_H_ */
