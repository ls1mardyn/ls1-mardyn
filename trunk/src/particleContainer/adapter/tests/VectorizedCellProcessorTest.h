/*
 * VectorizedCellProcessorTest.h
 *
 * @Date: 18.08.2014
 * @Author: tchipevn
 */

#ifndef VECTORIZEDCELLPROCESSORTEST_H_
#define VECTORIZEDCELLPROCESSORTEST_H_

#include "utils/TestWithSimulationSetup.h"

/**
 * This class tests the VectorizedCellProcessor, mostly against the Legacy one.
 *
 */
class VectorizedCellProcessorTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(VectorizedCellProcessorTest);

	TEST_METHOD(testForcePotentialCalculationU0);
	TEST_METHOD(testForcePotentialCalculationF0);

	TEST_METHOD(testLennardJonesVectorization);

	TEST_METHOD(testChargeChargeVectorization);
	TEST_METHOD(testChargeDipoleVectorization);
	TEST_METHOD(testChargeQuadrupoleVectorization);

	TEST_METHOD(testDipoleDipoleVectorization);
	TEST_METHOD(testDipoleQuadrupoleVectorization);

	TEST_METHOD(testQuadrupoleQuadrupoleVectorization);

	TEST_METHOD(testWaterVectorization);

	TEST_METHOD(testMultiComponentMultiPotentials);

	TEST_SUITE_END();

public:

	VectorizedCellProcessorTest();

	virtual ~VectorizedCellProcessorTest();

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

	/**
	 * Generic test routine for all electrostatic interactions.
	 * Which test is run is dependent on filename,
	 * i.e. what the scenario, which is being run, contains.
	 */
	void testElectrostaticVectorization(const char* filename, double ScenarioCutoff);

	/**
	 * Scenario, containing only charges.
	 */
	void testChargeChargeVectorization();

	/**
	 * Scenario, containing charges and dipoles.
	 */
	void testChargeDipoleVectorization();

	/**
	 * Scenario, containing charges and quadrupoles.
	 */
	void testChargeQuadrupoleVectorization();

	/**
	 * Scenario, containing only dipoles.
	 */
	void testDipoleDipoleVectorization();

	/**
	 * Scenario, containing dipoles and quadrupoles.
	 */
	void testDipoleQuadrupoleVectorization();

	/**
	 * Scenario, containing only quadrupoles.
	 */
	void testQuadrupoleQuadrupoleVectorization();

	/**
	 * Scenario, containing 1 LJ and 3 charges.
	 */
	void testWaterVectorization();

	/**
	 * Scenario, containing 2 components.
	 */
	void testMultiComponentMultiPotentials();

};
#endif /* VECTORIZEDCELLPROCESSORTEST_H_ */
