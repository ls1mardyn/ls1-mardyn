/*
 * MoleculeTest.h
 *
 * @Date: 04.02.2011
 * @Author: eckhardw
 */

#ifndef MOLECULETEST_H_
#define MOLECULETEST_H_

#include "utils/Testing.h"

class MoleculeTest : public utils::Test {

	TEST_SUITE(MoleculeTest);
	TEST_METHOD(testIsLessThan);
	TEST_SUITE_END();

public:

	MoleculeTest();

	virtual ~MoleculeTest();

	void testIsLessThan();

};

#endif /* MOLECULETEST_H_ */
