#pragma once

#include "../Testing.h"

class FixedSizeQueueTest : public utils::Test {
	TEST_SUITE(FixedSizeQueueTest);
	TEST_METHOD(testDefaultConstructed);
	TEST_METHOD(testFiveElements);
	TEST_METHOD(testZeroSize);
	TEST_METHOD(testGrowth);
	TEST_METHOD(testShrink);
	TEST_SUITE_END();

public:
	static void testDefaultConstructed();
	static void testFiveElements();
	static void testZeroSize();
	static void testGrowth();
	static void testShrink();
};
