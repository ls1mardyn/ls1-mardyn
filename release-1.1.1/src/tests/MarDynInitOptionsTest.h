/*
 * MarDynTest.h
 *
 * @Date: 12.07.2010
 * @Author: Wolfgang Eckhardt
 */

#ifndef MARDYNINITOPTIONSTEST_H_
#define MARDYNINITOPTIONSTEST_H_

#include "utils/Testing.h"

/**
 * This class tests, if the command-line options are all present and configured
 * the right way.
 */
class MarDynInitOptionsTest: public utils::Test {

	TEST_SUITE(MarDynInitOptionsTest);

	TEST_METHOD(testAllOptions);

	TEST_SUITE_END();

public:

	MarDynInitOptionsTest();

	virtual ~MarDynInitOptionsTest();

	void testAllOptions();
};

#endif /* MARDYNINITOPTIONSTEST_H_ */
