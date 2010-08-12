/*
 * CommonTest.cpp
 *
 * @Date: 21.05.2010
 * @Author: Wolfgang Eckhardt
 */

#include "CommonTest.h"
#include "Common.h"
#include "utils/Logger.h"

#include <iostream>

TEST_SUITE_REGISTRATION(CommonTest);

CommonTest::CommonTest() {
}

CommonTest::~CommonTest() {
}

void CommonTest::testGetTimeString() {
	std::string time = gettimestring();
	ASSERT_EQUAL((size_t)13, time.size());
	ASSERT_EQUAL(time[6], 'T');
}


void CommonTest::testAlignedNumber() {
//	ASSERT_EQUAL_MSG("One should be zero!", 1, 0);

	std::string result = aligned_number(123, 7, '.');
	std::string expected("....123");
	ASSERT_EQUAL_MSG("Align number 123 to 7 digits, filling char is .", expected, result);

	result = aligned_number(-2, 6, ' ');
	expected = "    -2";
	ASSERT_EQUAL_MSG("Align number -2 to 6 digits, filling char is ' '(space)", expected, result);
}
