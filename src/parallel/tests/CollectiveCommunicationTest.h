/*
 * CollectiveCommunicationTest.h
 *
 *  Created on: May 3, 2017
 *      Author: seckler
 */
#pragma once

#include "utils/Testing.h"

#include <vector>
class CollectiveCommunicationInterface;


class CollectiveCommunicationTest : public utils::Test {

	TEST_SUITE(CollectiveCommunicationTest);
	TEST_METHOD(testCollectiveCommBase);
	TEST_METHOD(testCollectiveCommunication);
	TEST_METHOD(testCollectiveCommunicationNonBlocking);
	TEST_SUITE_END();

public:

	CollectiveCommunicationTest();

	virtual ~CollectiveCommunicationTest();

	void testCollectiveCommBase();

	void testCollectiveCommunication();

	void testCollectiveCommunicationNonBlocking();

private:
	void testSingleIteration(CollectiveCommunicationInterface& collComm);
	int _rank;
	int _commSize;

};

