/**
 * \file
 * \brief
 * \author Johannes Heckl
 */

#include "AlignedArrayTest.h"
#include "../AlignedArray.h"

#include <cstdint>
#include <vector>

using namespace std;

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

void AlignedArrayTest::testAppending() {
	vector<double> a;
	for(int i = 0; i < 17; ++i) {
		a.push_back(static_cast<double>(i));
	}

	AlignedArray<double> aligned_a;
	for(int i = 0; i < 17; ++i) {
		aligned_a.appendValue(static_cast<double>(i), i);
	}

	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(&aligned_a[0]) % 64), 0);

	for(int i = 0; i < 17; ++i) {
		ASSERT_DOUBLES_EQUAL(a[i], aligned_a[i], 0.0);
	}
}

template<typename T>
void check(const vector<T>& a, const AlignedArray<T>& aligned_a, int i) {
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(&aligned_a[0]) % 64), 0);

	for (int j = 0; j < i; ++j) {
		mardyn_assert(a[j] == aligned_a[j]);
		ASSERT_DOUBLES_EQUAL(a[j], aligned_a[j], 0.0);
	}
}

void AlignedArrayTest::testIncreasingStorage() {
	vector<double> a;

	for(int i = 0; i < 111; ++i) {
		a.push_back(static_cast<double>(i));
	}

	AlignedArray<double> aligned_a;
	for(int i = 0; i < 11; ++i) {
		aligned_a.appendValue(static_cast<double>(i), i);
		check(a, aligned_a, i);
	}

	aligned_a.increaseStorage(11, 2);

	for(int i = 11; i < 22; ++i) {
		aligned_a.appendValue(static_cast<double>(i), i);
		check(a, aligned_a, i);
	}

	aligned_a.increaseStorage(22, 7);

	for(int i = 22; i < 111; ++i) {
		aligned_a.appendValue(static_cast<double>(i), i);
		check(a, aligned_a, i);
	}

	check(a, aligned_a, 111);
}
