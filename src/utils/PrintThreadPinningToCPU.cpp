/*
 * PrintThreadPinningToCPU.cpp
 *
 *  Created on: 26 Sep 2017
 *      Author: tchipevn
 */

#include "PrintThreadPinningToCPU.h"
#include "../WrapOpenMP.h"

#include <sched.h> /* int sched_getcpu(void); */
#include <iostream> /* intentionally writing to cout and not to global_log */

void PrintThreadPinningToCPU() {
	int mastercpu = sched_getcpu();

	std::cout << "Thread pinning:" << std::endl;
	std::cout << "Master thread running on " << mastercpu << std::endl;

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int cpu = sched_getcpu();
		const int myID = mardyn_get_thread_num();
		const int numThreads = mardyn_get_num_threads();

		for (int i = 0; i < numThreads; ++i) {
			if (i == myID) {
				std::cout << "	Thread with id " << myID << " is running on " << cpu << "." << std::endl;
			}

			#if defined(_OPENMP)
			#pragma omp barrier
			#endif
		}
	}
}
