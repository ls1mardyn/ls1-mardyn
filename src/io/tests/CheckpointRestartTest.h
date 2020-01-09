/*
 * CheckpointRestartTest.h
 *
 * Check whether a checkpoint can be successfully read again.
 *
 *  Created on: 11.08.2016
 *      Author: seckler
 */
#pragma once

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
