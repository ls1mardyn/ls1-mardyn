/*
 * MarDynTest.h
 *
 * @Date: 12.07.2010
 * @Author: Wolfgang Eckhardt
 */

#ifndef MARDYNINITOPTIONSTEST_H_
#define MARDYNINITOPTIONSTEST_H_

#include <cppunit/extensions/HelperMacros.h>

/**
 * This class tests, if the command-line options are all present and configured
 * the right way.
 */
class MarDynInitOptionsTest: public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(MarDynInitOptionsTest);

	CPPUNIT_TEST(testAllOptions);

	CPPUNIT_TEST_SUITE_END();

public:
	MarDynInitOptionsTest();
	virtual ~MarDynInitOptionsTest();

	void testAllOptions();

};

#endif /* MARDYNINITOPTIONSTEST_H_ */
