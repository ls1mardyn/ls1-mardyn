/*
 * WignerRotationTest.h
 *
 *  Created on: Jun 16, 2015
 *      Author: uwe
 */

#ifndef SRC_BHFMM_TESTS_WIGNERROTATIONTEST_H_
#define SRC_BHFMM_TESTS_WIGNERROTATIONTEST_H_

#include "utils/Testing.h"

/**
 * This class tests the VectorizedCellProcessor, mostly against the Legacy one.
 *
 */
class WignerRotationTest : public utils::Test {

	TEST_SUITE(WignerRotationTest);

	TEST_METHOD(testM2MWignerRotation);
	TEST_METHOD(testL2LWignerRotation);
	TEST_METHOD(testM2LWignerRotation);

	TEST_SUITE_END();

public:

	WignerRotationTest();

	virtual ~WignerRotationTest();

	/**
	 * Scenario, test M2M Operator based on Wigner rotation
	 */
	void testM2MWignerRotation();


	/**
	 * Scenario, test L2L Operator based on Wigner rotation
	 */
	void testL2LWignerRotation();

	/**
	 * Scenario, test M2L Operator based on Wigner rotation
	 */
	void testM2LWignerRotation();

};



#endif /* SRC_BHFMM_TESTS_WIGNERROTATIONTEST_H_ */
