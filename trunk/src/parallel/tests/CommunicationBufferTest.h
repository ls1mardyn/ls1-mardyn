/*
 * CommunicationBufferTest.h
 *
 *  Created on: 16 Oct 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARALLEL_TESTS_COMMUNICATIONBUFFERTEST_H_
#define SRC_PARALLEL_TESTS_COMMUNICATIONBUFFERTEST_H_

#include "utils/TestWithSimulationSetup.h"

class CommunicationBufferTest : public utils::TestWithSimulationSetup {
	TEST_SUITE(CommunicationBufferTest);
	TEST_METHOD(testEmplaceRead);
	TEST_METHOD(testHalo);
	TEST_METHOD(testLeaving);
	TEST_METHOD(testLeavingAndHalo);
	TEST_METHOD(testPackSendRecvUnpack);
	TEST_SUITE_END();

public:
	CommunicationBufferTest();

	virtual ~CommunicationBufferTest();

	void testEmplaceRead();

	void testHalo();

	void testLeaving();

	void testLeavingAndHalo();

	void testPackSendRecvUnpack();
};

#endif /* SRC_PARALLEL_TESTS_COMMUNICATIONBUFFERTEST_H_ */
