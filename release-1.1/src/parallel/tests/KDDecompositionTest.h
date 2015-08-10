/*
 * KDDecompositionTest.h
 *
 * @Date: 01.03.2012
 * @Author: eckhardw
 */

#ifndef KDDECOMPOSITIONTEST_H_
#define KDDECOMPOSITIONTEST_H_

#include "utils/TestWithSimulationSetup.h"
#include "parallel/KDNode.h"

class KDDecompositionTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(KDDecompositionTest);
	TEST_METHOD(testCompleteTreeInfo);
	TEST_SUITE_END();

public:

	KDDecompositionTest();

	virtual ~KDDecompositionTest();

	void testCompleteTreeInfo();

	/**
	 * Initial implementation of completeTreeInfo(). Kept to test the new / current
	 * implementation against it.
	 */
	void completeTreeInfo(KDNode*& root, KDNode*& ownArea, int ownRank);
};

#endif /* KDDECOMPOSITIONTEST_H_ */
