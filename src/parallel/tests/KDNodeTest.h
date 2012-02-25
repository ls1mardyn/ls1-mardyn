/*
 * KDNodeTest.h
 *
 *  Created on: Feb 24, 2012
 *      Author: eckhardw
 */

#ifndef KDNODETEST_H_
#define KDNODETEST_H_

#include "utils/Testing.h"

class KDNodeTest : public utils::Test {

	TEST_SUITE(KDNodeTest);
	TEST_METHOD(testEqual);
	TEST_METHOD(testBuildKDTree);
	TEST_METHOD(testFindAreaForProcess);
	TEST_SUITE_END();

public:

	KDNodeTest();

	virtual ~KDNodeTest();

	void testEqual();

	void testBuildKDTree();

	void testFindAreaForProcess();

};

#endif /* KDNODETEST_H_ */
