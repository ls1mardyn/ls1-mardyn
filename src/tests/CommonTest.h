/*
 * CommonTest.h
 *
 * @Date: 21.05.2010
 * @Author: Wolfgang Eckhardt
 */

#ifndef COMMONTEST_H_
#define COMMONTEST_H_

#include <cppunit/extensions/HelperMacros.h>

/**
 * Test functionality of module Common
 */
class CommonTest : public CppUnit::TestFixture{

	// declare testsuite commonTest
	CPPUNIT_TEST_SUITE(CommonTest);

	// add two methods which perform tests
	CPPUNIT_TEST(testGetTimeString);
	CPPUNIT_TEST(testAlignedNumber);

	// end suite declaration
	CPPUNIT_TEST_SUITE_END();


public:
	CommonTest();

	virtual ~CommonTest();

	void testGetTimeString();

	void testAlignedNumber();
};

#endif /* COMMONTEST_H_ */
