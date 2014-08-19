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

	// declare testsuite
	TEST_SUITE(AlignedArrayTest);
	// add methods which perform tests
	TEST_METHOD(testAlignment);
	// end suite declaration
	TEST_SUITE_END();

public:
	AlignedArrayTest();

	virtual ~AlignedArrayTest();

	void testAlignment();
};

#endif /* COMMONTEST_H_ */
