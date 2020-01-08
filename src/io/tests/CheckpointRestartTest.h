
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
	TEST_METHOD(testCheckpointRestartASCII);

	// add a method which perform test
	TEST_METHOD(testCheckpointRestartBinary);

	// end suite declaration
	TEST_SUITE_END();

public:
	CheckpointRestartTest() = default;
	virtual ~CheckpointRestartTest() = default;

	void testCheckpointRestartASCII();

	void testCheckpointRestartBinary();
private:

	void testCheckpointRestart(bool binary);
	unsigned long getGlobalParticleNumber(ParticleContainer* particleContainer);
};

#endif /* INPUTFILETEST_H_ */
