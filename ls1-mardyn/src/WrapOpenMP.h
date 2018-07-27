/*
 * WrapOpenMP.h
 *
 *  Created on: 15 Sep 2016
 *      Author: tchipevn
 */

#ifndef SRC_WRAPOPENMP_H_
#define SRC_WRAPOPENMP_H_

#include "utils/mardyn_assert.h"

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
	typedef omp_lock_t mardyn_lock_t;
	inline void mardyn_set_lock(mardyn_lock_t* l) {omp_set_lock(l);}
	inline void mardyn_init_lock(mardyn_lock_t* l) {omp_init_lock(l);}
	inline void mardyn_unset_lock(mardyn_lock_t* l) {omp_unset_lock(l);}
	inline void mardyn_destroy_lock(mardyn_lock_t* l) {omp_destroy_lock(l);}

#else
	inline int mardyn_get_thread_num()  {return 0;}
	inline int mardyn_get_num_threads() {return 1;}
	inline int mardyn_get_max_threads() {return 1;}
	typedef int mardyn_lock_t;
	inline void mardyn_set_lock(mardyn_lock_t* l) {
		mardyn_assert(*l == 0);
		*l = 1;
	}
	inline void mardyn_init_lock(mardyn_lock_t* l) {
		mardyn_assert(l != nullptr);
	}
	inline void mardyn_unset_lock(mardyn_lock_t* l) {
		mardyn_assert(*l == 1);
		*l = 0;
	}
	inline void mardyn_destroy_lock(mardyn_lock_t* l) {
		mardyn_assert(l != nullptr);
	}
#endif

#endif /* SRC_WRAPOPENMP_H_ */
