/*
 * CanonicalEnsembleTest.h
 *
 * @Date: 18.02.2011
 * @Author: eckhardw
 */

#ifndef CANONICALENSEMBLETEST_H_
#define CANONICALENSEMBLETEST_H_

#include "utils/TestWithSimulationSetup.h"

/**
 * Simple tests for the canonical ensemble.
 *
 * @todo In it's current implementation the Canonical Ensemble is only testable
 *       if the code is compiled without -DENABLE_MPI.
 *       For the parallel case the test case needed to set up the global_simulation somehow,
 *       as the ensemble needs it for the communication.
 */
class CanonicalEnsembleTest: public utils::TestWithSimulationSetup {

	TEST_SUITE(CanonicalEnsembleTest);
	TEST_METHOD(UpdateNumMoleculesSequential);
	TEST_METHOD(UpdateNumMoleculesParallel);
	TEST_SUITE_END();

public:

	CanonicalEnsembleTest();

	virtual ~CanonicalEnsembleTest();

	//! Test the update of the number of molecules
	void UpdateNumMoleculesSequential();

	void UpdateNumMoleculesParallel();

};

#endif /* CANONICALENSEMBLETEST_H_ */
