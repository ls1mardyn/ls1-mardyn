/*
 * PermutationTest.h
 *
 *  Created on: 2 Jun 2017
 *      Author: tchipevn
 */

#pragma once

#include "../Testing.h"

class RandomTest : public utils::Test {

	TEST_SUITE(RandomTest);
	TEST_METHOD(testRnd);
	TEST_SUITE_END();


public:
	RandomTest() {
		// TODO Auto-generated constructor stub

	}
	virtual ~RandomTest() {
		// TODO Auto-generated destructor stub
	}

	/**
	 *  I used this unit test to find out what Random::rnd() generates.
	 *  It seems to be a uniformly distributed random number in [0, 1].
	 *  And to check whether it can work in an OpenMP-parallel way.
	 *  Test does nothing except check that the numbers are in [0, 1].
	 *  */
	void testRnd();

	/**
	 * this method is not implemented!
	 */
//	void testGaussDeviate();
};
