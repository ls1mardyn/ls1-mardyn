/*
 * File:   NeighborAcquirerTest.h
 * Author: bierth, seckler
 *
 * Created on February 27, 2018, 5:01 PM
 */

#ifndef SRC_PARALLEL_TESTS_NEIGHBOURCOMMUNICATIONSCHEMETEST_H
#define SRC_PARALLEL_TESTS_NEIGHBOURCOMMUNICATIONSCHEMETEST_H

#include "utils/TestWithSimulationSetup.h"
#include "parallel/ZonalMethods/FullShell.h"
#include "parallel/NeighbourCommunicationScheme.h"


class NeighborAcquirerTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(NeighborAcquirerTest);
	TEST_METHOD(testShiftIfNecessary);
	TEST_METHOD(testOverlap);
	TEST_METHOD(testIOwnThis);
	TEST_METHOD(testCorrectNeighborAcquisition);
	TEST_SUITE_END();

	public:
		NeighborAcquirerTest();
		~NeighborAcquirerTest();
		void testShiftIfNecessary();
		void testOverlap();
		void testIOwnThis();
		void testCorrectNeighborAcquisition();
	private:
		FullShell* _fullShell;
};


#endif /* SRC_PARALLEL_TESTS_NEIGHBOURCOMMUNICATIONSCHEMETEST_H */

