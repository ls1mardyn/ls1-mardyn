/*
 * PermutationTest.cpp
 *
 *  Created on: 2 Jun 2017
 *      Author: tchipevn
 */

#include "RandomTest.h"
#include "../Random.h"
#include "../../WrapOpenMP.h"

#include <iostream>

TEST_SUITE_REGISTRATION(RandomTest);

void RandomTest::testRnd() {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		int myID = mardyn_get_thread_num();
		int seed = myID;
		Random threadPrivateRNG = Random(seed);

		int numSamples = 1000;
//		const int numBins = 10000;
//		std::array<int, numBins> bins;
		for (int i = 0; i < numSamples; ++i) {

			float sample = threadPrivateRNG.rnd();

			ASSERT_TRUE(sample <= 1.0f);
			ASSERT_TRUE(sample >= 0.0f);

			// find correct bin
			int bin = int(sample * 10);
			ASSERT_TRUE(bin >= 0);
			ASSERT_TRUE(bin < 10);
//
//#pragma omp critical
//			std::cout << myID << " " << threadPrivateRNG.rnd() << std::endl;
		}
	}
}
