/**
 * \file
 * \brief
 * \author Johannes Heckl
 */

#include "AlignedArrayTest.h"
#include "../AlignedArray.h"
#include <stdint.h>

TEST_SUITE_REGISTRATION(AlignedArrayTest);

AlignedArrayTest::AlignedArrayTest() {
}

AlignedArrayTest::~AlignedArrayTest() {
}

void AlignedArrayTest::testAlignment() {
	Log::global_log->info() << "Testing AlignedArray." << std::endl ;

	const size_t length = 16;
	AlignedArray<double, 16> a(length);
	AlignedArray<double, 32> b(length);
	AlignedArray<double, 64> c(length);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(&a[0]) % 16), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(&b[0]) % 32), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(&c[0]) % 64), 0);
}
