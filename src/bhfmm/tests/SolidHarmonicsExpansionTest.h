/*
 * SolidHarmonicsExpansionTest.h
 *
 *  Created on: 2 Feb 2018
 *      Author: tchipevn
 */

#ifndef SRC_BHFMM_TESTS_SOLIDHARMONICSEXPANSIONTEST_H_
#define SRC_BHFMM_TESTS_SOLIDHARMONICSEXPANSIONTEST_H_

#include "utils/Testing.h"

class SolidHarmonicsExpansionTest : public utils::Test {

	TEST_SUITE(SolidHarmonicsExpansionTest);
	TEST_METHOD(test_P2M_M2P);
	TEST_METHOD(test_P2L_L2P);
	TEST_SUITE_END();

public:
	SolidHarmonicsExpansionTest();

	virtual ~SolidHarmonicsExpansionTest();

	void test_P2M_M2P();

	void test_P2L_L2P();
};

#endif /* SRC_BHFMM_TESTS_SOLIDHARMONICSEXPANSIONTEST_H_ */
