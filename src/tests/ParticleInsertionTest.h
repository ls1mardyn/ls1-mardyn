/*
 * ParticleInsertionTest.h
 * Tests usher algorithm
 *
 *  Created on: Jun 18, 2012
 *      Author: tijana
 */

#ifndef PARTICLEINSERTIONTEST_H_
#define PARTICLEINSERTIONTEST_H_

#include "utils/TestWithSimulationSetup.h"
#include "Simulation.h"
#include "molecules/Molecule.h"
#include "particleContainer/LinkedCells.h"
#include "ParticleInsertion.h"
#include <sstream>

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

class ParticleInsertionTest: public utils::TestWithSimulationSetup {
	TEST_SUITE( ParticleInsertionTest);
	//TEST_METHOD( testRotation);
	//TEST_METHOD( testTranslationAndRotation);
	//TEST_METHOD( testParameterSetup);
	//TEST_METHOD( testParametersFullStudy);
	TEST_SUITE_END();
public:
	ParticleInsertionTest();
	virtual ~ParticleInsertionTest();

	void testRotation();

	void testTranslationAndRotation();

	void testParameterSetup();

	void testParametersFullStudy();

	void readParamFile(string file_name, int* maxIter, int* maxRestarts, int* maxRotations, double* tolerance,
			double* maxAngle, double* maxAllowedAngle, double* minAngle, bool* largeStepsizeOnOverlap,
			bool* restartIfIncreases);
};

#endif /* PARTICLEINSERTIONTEST_H_ */
