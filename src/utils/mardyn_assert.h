/*
 * mardyn_assert.h
 *
 *  Created on: 8 Feb 2017
 *      Author: tchipevn
 */

#ifndef SRC_UTILS_MARDYN_ASSERT_H_
#define SRC_UTILS_MARDYN_ASSERT_H_

#include <string>
#include <sstream>

#include "Logger.h"

// Macro to wrap mardyn_exit and pass the caller file and line
#define MARDYN_EXIT(exit_message) mardyn_exit(exit_message, __FILE__, __LINE__)

inline void mardyn_exit(const std::ostringstream & exit_message,
						const char* file, const int line) {
	Log::global_log->error_always_output()
		<< "Exit called in file `" << file << ":" << line << "`" << std::endl;
	Log::global_log->error_always_output()
		<< exit_message << std::endl;
#ifdef ENABLE_MPI
	// terminate all mpi processes and return exitcode
	MPI_Abort(MPI_COMM_WORLD, code);
#else
	// call global abort - this stops the debugger at the right spot.
	::abort();
#endif
}

inline void __mardyn_assert__(const char * expr, const char* file, int line) {
	std::ostringstream error_message;
	error_message << "Assertion \"" << expr << "\" failed at " << file << ":" << line << std::endl;
	MARDYN_EXIT(error_message);
}

#ifdef NDEBUG
#define mardyn_assert(EXPRESSION) ((void)0)
#else
#define mardyn_assert(EXPRESSION) ((EXPRESSION) ? (void)0 : __mardyn_assert__(#EXPRESSION, __FILE__, __LINE__))
#endif

#endif /* SRC_UTILS_MARDYN_ASSERT_H_ */
