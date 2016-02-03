
/*
 * inputFileTest.h
 *
 *  Created on: 01.05.2012
 *      Author: yutaka
 */
#ifdef SUPPORT_GENERATOR
#ifndef INPUTFILETEST_H_
#define INPUTFILETEST_H_

#include "utils/TestWithSimulationSetup.h"

class InputFileTest:public utils::TestWithSimulationSetup {

	// declare testsuite inputFileTest
	TEST_SUITE(InputFileTest);

	// add a method which perform test
	TEST_METHOD(testRemoveMomentum);

	// end suite declaration
	TEST_SUITE_END();

public:
	InputFileTest();
	virtual ~InputFileTest();
	/*
	 * testRemoveMomentum tests if removeMomentum in MDGenerator works properly or not.
	 */
	void testRemoveMomentum();
};

#endif /* INPUTFILETEST_H_ */

#endif/*SUPPORT_GENERATOR*/
