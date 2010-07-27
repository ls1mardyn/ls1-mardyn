/*
 * CommonTest.cpp
 *
 * @Date: 21.05.2010
 * @Author: eckhardw
 */

#include "CommonTest.h"
#include "Common.h"
#include "utils/Logger.h"

CPPUNIT_TEST_SUITE_REGISTRATION(CommonTest);

CommonTest::CommonTest() {
}

CommonTest::~CommonTest() {
}

void CommonTest::testGetTimeString() {
	//Log::Logger log(Log::All);
	//log.debug() << "Testing getTimeString" << std::endl;
	std::string time = gettimestring();
	CPPUNIT_ASSERT_EQUAL((size_t)13, time.size());
	CPPUNIT_ASSERT_EQUAL(time[6], 'T');
}


void CommonTest::testAlignedNumber() {
	//CPPUNIT_ASSERT_EQUAL_MESSAGE("One should be zero!", 1, 0);

	std::string result = aligned_number(123, 7, '.');
	std::string expected("....123");
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Align number 123 to 7 digits, filling char is .", expected, result);

	result = aligned_number(-2, 6, ' ');
	expected = "    -2";
	CPPUNIT_ASSERT_EQUAL_MESSAGE("Align number -2 to 6 digits, filling char is ' '(space)", expected, result);
}
