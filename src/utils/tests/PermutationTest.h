/*
 * PermutationTest.h
 *
 *  Created on: 2 Jun 2017
 *      Author: tchipevn
 */

#pragma once

#include "../Testing.h"

class PermutationTest : public utils::Test {

	TEST_SUITE(PermutationTest);
	TEST_METHOD(testPermutations);
	TEST_SUITE_END();


public:
	PermutationTest() {
		// TODO Auto-generated constructor stub

	}
	virtual ~PermutationTest() {
		// TODO Auto-generated destructor stub
	}

	void testPermutations();
};
