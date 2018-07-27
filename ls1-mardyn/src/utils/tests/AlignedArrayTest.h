/**
 * \file
 * \brief Tests for AlignedArray.
 * \author Johannes Heckl
 */

#ifndef ALIGNEDARRAYTEST_H_
#define ALIGNEDARRAYTEST_H_

#include "../Testing.h"

/**
 * \brief Tests for AlignedArray.
 * \author Johannes Heckl
 */
class AlignedArrayTest : public utils::Test {

	TEST_SUITE(AlignedArrayTest);
	TEST_METHOD(testAlignment);
	TEST_METHOD(testAppending);
	TEST_METHOD(testIncreasingStorage);
	TEST_SUITE_END();

public:
	AlignedArrayTest();

	virtual ~AlignedArrayTest();

	void testAlignment();

	void testAppending();

	void testIncreasingStorage();
};

#endif /* COMMONTEST_H_ */
