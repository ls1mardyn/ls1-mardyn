/*
 * AlignedArrayTripletTest.cpp
 *
 *  Created on: 20 Jun 2017
 *      Author: tchipevn
 */

#include "AlignedArrayTripletTest.h"
#include "../AlignedArrayTriplet.h"

#include <cstdint>
#include <vector>


TEST_SUITE_REGISTRATION(AlignedArrayTripletTest);

AlignedArrayTripletTest::AlignedArrayTripletTest() {
}

AlignedArrayTripletTest::~AlignedArrayTripletTest() {
}


void AlignedArrayTripletTest::testAlignment() {
	Log::global_log->info() << "Testing AlignedArrayTriplet." << std::endl ;

	const size_t length = 16;
	AlignedArrayTriplet<double> a(length);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(a.xBegin()) % 64), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(a.yBegin()) % 64), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(a.zBegin()) % 64), 0);
}

void AlignedArrayTripletTest::testAppending() {
	std::vector<double> x, y, z;
	for(int i = 0; i < 17; ++i) {
		x.push_back(static_cast<double>(i));
		y.push_back(static_cast<double>(i+1));
		z.push_back(static_cast<double>(i+2));
	}

	AlignedArrayTriplet<double> A;
	for(int i = 0; i < 17; ++i) {
		A.appendValueTriplet(static_cast<double>(i), static_cast<double>(i+1), static_cast<double>(i+2), i);
	}

	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(A.xBegin()) % 64), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(A.yBegin()) % 64), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(A.zBegin()) % 64), 0);

	for(int i = 0; i < 17; ++i) {
		ASSERT_DOUBLES_EQUAL(x[i], A.x(i), 0.0);
		ASSERT_DOUBLES_EQUAL(y[i], A.y(i), 0.0);
		ASSERT_DOUBLES_EQUAL(z[i], A.z(i), 0.0);
	}
}

template<typename T>
void check(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& z, const AlignedArrayTriplet<T>& A, int i) {
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(A.xBegin()) % 64), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(A.yBegin()) % 64), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(A.zBegin()) % 64), 0);

	for (int j = 0; j < i; ++j) {
		mardyn_assert(x[j] == A.x(j));
		mardyn_assert(y[j] == A.y(j));
		mardyn_assert(z[j] == A.z(j));

		ASSERT_DOUBLES_EQUAL(x[j], A.x(j), 0.0);
		ASSERT_DOUBLES_EQUAL(y[j], A.y(j), 0.0);
		ASSERT_DOUBLES_EQUAL(z[j], A.z(j), 0.0);
	}
}

void AlignedArrayTripletTest::testIncreasingStorage() {
	std::vector<double> x, y, z;

	for(int i = 0; i < 111; ++i) {
		x.push_back(static_cast<double>(i));
		y.push_back(static_cast<double>(-i-1));
		z.push_back(static_cast<double>(i+2));
	}

	AlignedArrayTriplet<double> A;
	for(int i = 0; i < 11; ++i) {
		A.appendValueTriplet(static_cast<double>(i), static_cast<double>(-i-1), static_cast<double>(i+2), i);
		check(x, y, z, A, i);
	}

	A.increaseStorage(11, 2);

	for(int i = 11; i < 22; ++i) {
		A.appendValueTriplet(static_cast<double>(i), static_cast<double>(-i-1), static_cast<double>(i+2), i);
		check(x, y, z, A, i);
	}

	A.increaseStorage(22, 7);

	for(int i = 22; i < 111; ++i) {
		A.appendValueTriplet(static_cast<double>(i), static_cast<double>(-i-1), static_cast<double>(i+2), i);
		check(x, y, z, A, i);
	}

	check(x, y, z, A, 111);
}
