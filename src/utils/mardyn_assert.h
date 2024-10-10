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
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "Logger.h"

// Macro to wrap mardyn_exit and pass the caller file and line
#define MARDYN_EXIT(exit_message) mardyn_exit(exit_message, __func__, __FILE__, __LINE__)

inline void mardyn_exit(const std::string & exit_message,
						const char* function, const char* filepath, const int line) {
	
	// Only print the file path relative to the "/src/" directory
	// The following code extracts this relative path from "filepath"
	// Search until last "/src/" was found to avoid user-specific "src" dirs
	const char* filepath_truncated = nullptr;
	const char* temp = filepath;
	while ((temp = std::strstr(temp, "/src/")) != nullptr) {
		filepath_truncated = temp;
		temp++;
	}
	if (filepath_truncated == nullptr) {
		// Print absolute path if "/src/" not found
		filepath_truncated = filepath;
	} else {
  		filepath_truncated++;  // +1 to ignore the first "/"
	}

	// Print code location from which MARDYN_EXIT() was called
	Log::global_log->error()
		<< "Exit called from function `" << function << "`"
		<< " in file `" << filepath_truncated << ":" << line
		<< "` with message:" << std::endl;

	// Print exit message line by line to always have Logger output
	std::stringstream ss(exit_message);
	std::string exit_message_line;
	// Split exit message by "\n" and print via Logger
	while (std::getline(ss, exit_message_line, '\n')) {
		Log::global_log->error() << exit_message_line << std::endl;
	}

#ifdef ENABLE_MPI
	// terminate all mpi processes and return exitcode
	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#else
	// call global abort - this stops the debugger at the right spot.
	::abort();
#endif
}

inline void __mardyn_assert__(const char * expr, const char* file, int line) {
	std::ostringstream error_message;
	error_message << "Assertion \"" << expr << "\" failed at " << file << ":" << line << std::endl;
	MARDYN_EXIT(error_message.str());
}

#ifdef NDEBUG
#define mardyn_assert(EXPRESSION) ((void)0)
#else
#define mardyn_assert(EXPRESSION) ((EXPRESSION) ? (void)0 : __mardyn_assert__(#EXPRESSION, __FILE__, __LINE__))
#endif

#endif /* SRC_UTILS_MARDYN_ASSERT_H_ */
