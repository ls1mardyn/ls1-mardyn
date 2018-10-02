/**
 * @file ZonalMethodTest.h
 * @author seckler
 * @date 02.10.18
 */

#pragma once

#include <utils/Testing.h>

class ZonalMethodTest : public utils::Test {
	TEST_SUITE(ZonalMethodTest);
	TEST_METHOD(testES);
	TEST_SUITE_END();
public:
	ZonalMethodTest() = default;
	~ZonalMethodTest() override = default;
	void testES();
};




