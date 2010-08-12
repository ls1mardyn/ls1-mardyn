/*
 * CommonTest.h
 *
 * @Date: 21.05.2010
 * @Author: Wolfgang Eckhardt
 */

#ifndef COMMONTEST_H_
#define COMMONTEST_H_

#include "utils/Testing.h"

/**
 * Test functionality of module Common
 */
class CommonTest : public utils::Test {

	// declare testsuite commonTest
	TEST_SUITE(CommonTest);

	// add two methods which perform tests
	TEST_METHOD(testGetTimeString);
	TEST_METHOD(testAlignedNumber);

	// end suite declaration
	TEST_SUITE_END();

public:
	CommonTest();

	virtual ~CommonTest();

	void testGetTimeString();

	void testAlignedNumber();
};

#endif /* COMMONTEST_H_ */
