/*
 * VCP1CLJWRTest.h
 *
 *  Created on: 13 Apr 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_TESTS_VCP1CLJWRTEST_H_
#define SRC_PARTICLECONTAINER_ADAPTER_TESTS_VCP1CLJWRTEST_H_

#include "utils/TestWithSimulationSetup.h"

class VCP1CLJWRTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(VCP1CLJWRTest);

	TEST_METHOD(testForcePotentialCalculationU0);
	TEST_METHOD(testForcePotentialCalculationF0);
	TEST_METHOD(testLennardJonesVectorization);

	TEST_SUITE_END();


public:
	VCP1CLJWRTest();
	virtual ~VCP1CLJWRTest();

	/**
	 * p3 --- p4
	 *  |      |
	 *  |      |
	 * p1 --- p2
	 *
	 * The particles are arranged on a rectangle so that the pairs p1-p2, p1-p3,
	 * p2-p4 and p3-p4 interact. Values are chosen so that potential is zero, force is
	 * exactly 1.0. (epsilon = 1.0, sigma = 1.0, r_{ab} = 1.0)
	 *
	 * See tests/integrationTests/ForceCalculation
	 */
	void testForcePotentialCalculationU0();

	/**
	 * Same setup as above, however potential should be U= and force F=0.
	 *
	 * @see testForcePotentialCalculationU0
	 */
	void testForcePotentialCalculationF0();

	/**
	 * Run first the Legacy- and then the VectorizedCellProcessor on the same
	 * input file. Then verify all forces and torques on each molecule, as well as
	 * the potential and virial.
	 * Input file contains only multi-centered Lennard-Jones components.
	 */
	void testLennardJonesVectorization();

};

#endif /* SRC_PARTICLECONTAINER_ADAPTER_TESTS_VCP1CLJWRTEST_H_ */
