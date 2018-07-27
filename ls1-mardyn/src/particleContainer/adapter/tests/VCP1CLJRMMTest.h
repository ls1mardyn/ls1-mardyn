#ifndef SRC_PARTICLECONTAINER_ADAPTER_TESTS_VCP1CLJRMMTEST_H_
#define SRC_PARTICLECONTAINER_ADAPTER_TESTS_VCP1CLJRMMTEST_H_

#include "utils/TestWithSimulationSetup.h"

class VCP1CLJRMMTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(VCP1CLJRMMTest);

	TEST_METHOD(testForcePotentialCalculationU0);
	TEST_METHOD(testForcePotentialCalculationF0);
	TEST_METHOD(testProcessCell);
	TEST_METHOD(testProcessCellPair);

	TEST_SUITE_END();


public:
	VCP1CLJRMMTest();
	virtual ~VCP1CLJRMMTest();

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
	 * input file. Then verify all forces and the potential computed on the first cell.
	 * Input file contains only single-centered Lennard-Jones components.
	 */
	void testProcessCell();

	/**
	 * Run first the Legacy- and then the VectorizedCellProcessor on the same
	 * input file. Then verify all forces and the potential computed on the first cell.
	 * Input file contains only single-centered Lennard-Jones components.
	 */
	void testProcessCellPair();


};

#endif /* SRC_PARTICLECONTAINER_ADAPTER_TESTS_VCP1CLJRMMTEST_H_ */
