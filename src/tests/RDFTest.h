/*
 * RDFTest.h
 *
 * @Date: 15.02.2011
 * @Author: eckhardw
 */

#ifndef RDFTEST_H_
#define RDFTEST_H_

#include "utils/TestWithSimulationSetup.h"


class RDFTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(RDFTest);
	TEST_METHOD(testRDFCountSequential12_LinkedCell);
	TEST_METHOD(testRDFCountSequential12_AdaptiveCell);
	TEST_METHOD(testRDFCountLinkedCell);
	TEST_METHOD(testRDFCountAdaptiveCell);
	TEST_SUITE_END();

public:

	RDFTest();

	virtual ~RDFTest();

	/**
	 * test the RDF for a system of 2x2x3 particles, for the sequential implementation,
	 * so it's not decomposed into several subdomains.
	 *
	 * First determine the particle pairs with the halo being empty . All particles are
	 * located in a regular cartesian grid with mesh width 1.0, the first
	 * particle being located at (1.0/1.0/1.0).
	 *
	 * As we use a cutoff-radius of 1.8 with 100 bins, this implys that in the grid
	 * occur 3 distances for particle pairs: 1.0, sqrt(2.0), and sqrt(3.0), so that
	 * only particle pairs in bins 55 (count = 20), 78 (count = 22) and 96 (count = 8)
	 * can be counted.
	 *
	 * Then we fill the halo and determine the RDF again. Due to the layout of the
	 * domain, only particles with x3==3.0 are copied into the halo with x3==-0.5,
	 * so for them the only distance to particles within the cutoff is 1.5, corresponding
	 * to bin 83 (count = 4).
	 */
	void testRDFCountSequential12_LinkedCell();

	void testRDFCountSequential12_AdaptiveCell();


	/**
	 * Performs the actual test, parameterized by the particleContainer. Thus we
	 * make sure that the RDF is determined valid for different implementations
	 * of particle containers.
	 *
	 * @todo Actually the test of different particle containers does not belong here!
	 *       But as we don't have any tests for particle containers, I guess this is ok.
	 */
	void testRDFCountSequential12(ParticleContainer* moleculeContainer);

	/**
	 * Test with 12x12x12 molecules, also parallel. In this case always the halo
	 * has to be populated (otherwise the number of subdomains would be wrong as we
	 * wouldn't calculate number of pairs over domain boarders.
	 */
	void testRDFCountLinkedCell();

	void testRDFCountAdaptiveCell();

	void testRDFCount(ParticleContainer* container);
};

#endif /* RDFTEST_H_ */
