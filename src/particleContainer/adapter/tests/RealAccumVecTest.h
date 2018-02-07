/*
 * RealAccumVecTest.h
 *
 *  Created on: 10 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_TESTS_REALACCUMVECTEST_H_
#define SRC_PARTICLECONTAINER_ADAPTER_TESTS_REALACCUMVECTEST_H_

#include "utils/Testing.h"

class RealAccumVecTest : public utils::Test  {
	TEST_SUITE(RealAccumVecTest);

	TEST_METHOD(testConvert);

	TEST_SUITE_END();

public:
	RealAccumVecTest();
	virtual ~RealAccumVecTest();

	void setUp();
	void tearDown();

	void testConvert();

};

#endif /* SRC_PARTICLECONTAINER_ADAPTER_TESTS_REALACCUMVECTEST_H_ */
