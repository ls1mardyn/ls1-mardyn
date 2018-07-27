/*
 * ConcatenatedAlignedArrayRMMTest.h
 *
 *  Created on: 30 Sep 2017
 *      Author: tchipevn
 */

#ifndef SRC_UTILS_TESTS_CONCATENATEDALIGNEDARRAYRMMTEST_H_
#define SRC_UTILS_TESTS_CONCATENATEDALIGNEDARRAYRMMTEST_H_

#include "../Testing.h"

class ConcatenatedAlignedArrayRMMTest: public utils::Test {
	TEST_SUITE(ConcatenatedAlignedArrayRMMTest);
	TEST_METHOD(testAlignment);
	TEST_METHOD(testGetters);
	TEST_METHOD(testZero);
	TEST_METHOD(testAppending);
	TEST_METHOD(testIncreasingStorage);
	TEST_SUITE_END();

public:
	ConcatenatedAlignedArrayRMMTest();

	virtual ~ConcatenatedAlignedArrayRMMTest();

	void testAlignment();

	void testGetters();

	void testZero();

	void testAppending();

	void testIncreasingStorage();
};

#endif /* SRC_UTILS_TESTS_CONCATENATEDALIGNEDARRAYRMMTEST_H_ */
