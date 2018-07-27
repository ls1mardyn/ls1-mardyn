/*
 * AlignedArrayTripletTest.h
 *
 *  Created on: 20 Jun 2017
 *      Author: tchipevn
 */

#ifndef SRC_UTILS_TESTS_ALIGNEDARRAYTRIPLETTEST_H_
#define SRC_UTILS_TESTS_ALIGNEDARRAYTRIPLETTEST_H_

#include "../Testing.h"

class AlignedArrayTripletTest: public utils::Test {
	TEST_SUITE(AlignedArrayTripletTest);
	TEST_METHOD(testAlignment);
	TEST_METHOD(testAppending);
	TEST_METHOD(testIncreasingStorage);
	TEST_SUITE_END();

public:
	AlignedArrayTripletTest();

	virtual ~AlignedArrayTripletTest();

	void testAlignment();

	void testAppending();

	void testIncreasingStorage();
};

#endif /* SRC_UTILS_TESTS_ALIGNEDARRAYTRIPLETTEST_H_ */
