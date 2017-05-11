
/*
 * inputFileTest.h
 *
 *  Created on: 01.05.2012
 *      Author: yutaka
 */
#ifndef INPUTFILETEST_H_
#define INPUTFILETEST_H_

#include "utils/TestWithSimulationSetup.h"

class CheckpointRestartTest : public utils::TestWithSimulationSetup {

	// declare testsuite inputFileTest
	TEST_SUITE(CheckpointRestartTest);

	// add a method which perform test
	TEST_METHOD(testCheckpointRestart);

	// end suite declaration
	TEST_SUITE_END();

public:
	CheckpointRestartTest();
	virtual ~CheckpointRestartTest();
	/*
	 * testRemoveMomentum tests if removeMomentum in MDGenerator works properly or not.
	 */
	void testCheckpointRestart();
};

#endif /* INPUTFILETEST_H_ */
