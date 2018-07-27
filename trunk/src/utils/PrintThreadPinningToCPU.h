/*
 * PrintThreadPinningToCPU.h
 *
 *  Created on: 26 Sep 2017
 *      Author: tchipevn
 */

#ifndef SRC_UTILS_PRINTTHREADPINNINGTOCPU_H_
#define SRC_UTILS_PRINTTHREADPINNINGTOCPU_H_

/**
 * this reproduces the functionality of Intel's KMP_AFFINITY=verbose
 * which is unfortunately unavailable in GNU's OpenMP implementation.
 */
void PrintThreadPinningToCPU();

#endif /* SRC_UTILS_PRINTTHREADPINNINGTOCPU_H_ */
