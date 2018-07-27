/*
 *  Created on: 30 May 2017
 *      Author: Micha Mueller
 */

#ifndef REALCALCVECTEST_H
#define REALCALCVECTEST_H

#include "utils/Testing.h"

class RealCalcVecTest : public utils::Test {

	TEST_SUITE(RealCalcVecTest);

	TEST_METHOD(testFastReciprocalMask);
	TEST_METHOD(testFastReciprocSqrtMask);

	TEST_SUITE_END();

public:
	RealCalcVecTest();
	virtual ~RealCalcVecTest();

	void setUp();

	void tearDown();

	/**
	 * Test if the results from fastReciprocal_mask() are the same as if the primitive approach was used
	 */
	void testFastReciprocalMask();

	/**
	 * Test if the results from fastReciprocSqrt_mask() are the same as if the primitive approach was used
	 */
	void testFastReciprocSqrtMask();

};

#endif /* REALCALCVECTEST_H */
