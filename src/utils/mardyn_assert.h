/*
 * mardyn_assert.h
 *
 *  Created on: 8 Feb 2017
 *      Author: tchipevn
 */

#ifndef SRC_UTILS_MARDYN_ASSERT_H_
#define SRC_UTILS_MARDYN_ASSERT_H_

#include "Logger.h"

inline void mardyn_exit(int code) {
#ifdef ENABLE_MPI
	// terminate all mpi processes and return exitcode
	MPI_Abort(MPI_COMM_WORLD, code);
#else
	// call global exit
	::exit(code);
#endif
}

inline void __mardyn_assert__(const char * expr, const char* file, int line) {
	Log::global_log->error_always_output() << "Assertion \"" << expr << "\" failed in file " << file << " on line " << line << std::endl;
	mardyn_exit(1);
}

#ifdef NDEBUG
#define mardyn_assert(EXPRESSION) ((void)0)
#else
#define mardyn_assert(EXPRESSION) ((EXPRESSION) ? (void)0 : __mardyn_assert__(#EXPRESSION, __FILE__, __LINE__))
#endif

#endif /* SRC_UTILS_MARDYN_ASSERT_H_ */
