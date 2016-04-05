/*
 * UnorderedVectorTest.h
 *
 *  Created on: Nov 19, 2015
 *      Author: tchipevn
 */

#ifndef UNORDEREDVECTORTEST_H_
#define UNORDEREDVECTORTEST_H_

#include "../Testing.h"

/**
 * \brief Test correctness of UnorderedVector fast removal.
 * \author Nikola Tchipev
 */
class UnorderedVectorTest: public utils::Test {
	TEST_SUITE(UnorderedVectorTest);
	TEST_METHOD(testFastRemovalInt);
	TEST_METHOD(testFastRemovalMoleculePointer);
	TEST_SUITE_END();

public:
	UnorderedVectorTest() {}
	virtual ~UnorderedVectorTest() {}

	void testFastRemovalInt();
	void testFastRemovalMoleculePointer();
};

#endif /* UNORDEREDVECTORTEST_H_ */
