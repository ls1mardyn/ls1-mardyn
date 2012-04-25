/*
 * ForceCalculationTest.h
 *
 * @Date: 22.04.2012
 * @Author: eckhardw
 */

#ifndef FORCECALCULATIONTEST_H_
#define FORCECALCULATIONTEST_H_

#include "utils/TestWithSimulationSetup.h"

/**
 * This class tests the iteration over the 4 1CLJ-molecules and the calculation
 * of the force and potential, respectively.
 */
class ForceCalculationTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(ForceCalculationTest);
	TEST_METHOD(testForcePotentialCalculation);
	TEST_SUITE_END();

public:

	ForceCalculationTest();

	virtual ~ForceCalculationTest();

	/**
	 * p3 --- p4
	 *  |      |
	 *  |      |
	 * p1 --- p2
	 *
	 * The particles are arranged on a rectangle so that the pairs p1-p2, p1-p3,
	 * p2-p4 and p3-p4 interact. Values are chosen so that potential and force are zero.
	 * (epsilon = 1.0, sigma = 1.0, r_{ab} = 1.0)
	 */
	void testForcePotentialCalculation();
};

#endif /* FORCECALCULATIONTEST_H_ */
