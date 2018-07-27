/*
 * MemoryProfilerTest.h
 *
 *  Created on: May 9, 2017
 *      Author: seckler
 */
#pragma once

#include "utils/Testing.h"

class MemoryProfilerTest: public utils::Test {

	// declare testsuite
TEST_SUITE(MemoryProfilerTest);
	// add methods which perform tests
	TEST_METHOD(testMemProfiler);
	// end suite declaration
	TEST_SUITE_END
	();

public:
	MemoryProfilerTest() {
	}
	virtual ~MemoryProfilerTest() {
	}
	void testMemProfiler();
};

