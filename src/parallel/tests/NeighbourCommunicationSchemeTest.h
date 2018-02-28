/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   NeighbourCommunicationSchemeTest.h
 * Author: bierth
 *
 * Created on February 27, 2018, 5:01 PM
 */

#ifndef SRC_PARALLEL_TESTS_NEIGHBOURCOMMUNICATIONSCHEMETEST_H
#define SRC_PARALLEL_TESTS_NEIGHBOURCOMMUNICATIONSCHEMETEST_H

#include "utils/TestWithSimulationSetup.h"
#include "parallel/ZonalMethods/FullShell.h"
#include "parallel/NeighbourCommunicationScheme.h"

class NeighbourCommunicationSchemeTest : public utils::TestWithSimulationSetup {
	
	TEST_SUITE(NeighbourCommunicationSchemeTest);
	TEST_METHOD(testShiftIfNecessary);
	TEST_METHOD(testOverlap);
	TEST_METHOD(testIOwnThis);
	TEST_SUITE_END();
	
	public:
		NeighbourCommunicationSchemeTest();
		~NeighbourCommunicationSchemeTest();
		void testShiftIfNecessary();
		void testOverlap();
		void testIOwnThis();
		
	private:
		FullShell* _fullShell;
		DirectNeighbourCommunicationScheme* _directScheme;
	
};


#endif /* SRC_PARALLEL_TESTS_NEIGHBOURCOMMUNICATIONSCHEMETEST_H */

