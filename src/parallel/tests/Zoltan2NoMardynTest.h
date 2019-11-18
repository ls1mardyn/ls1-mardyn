/**
 * @file Zoltan2NoMardynTest.h
 * @author seckler
 * @date 18.11.19
 */

#pragma once

#include <utils/Testing.h>

class Zoltan2NoMardynTest : public utils::Test {
	TEST_SUITE(Zoltan2NoMardynTest);
	TEST_METHOD(multiJaggedTest);
	TEST_SUITE_END();
public:
	Zoltan2NoMardynTest() = default;
	~Zoltan2NoMardynTest() override = default;

	void multiJaggedTest();
};




