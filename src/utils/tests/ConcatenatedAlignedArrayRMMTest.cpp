/*
 * ConcatenatedAlignedArrayRMMTest.cpp
 *
 *  Created on: 20 Jun 2017
 *      Author: tchipevn
 */

#include "ConcatenatedAlignedArrayRMMTest.h"
#include "../ConcatenatedAlignedArrayRMM.h"

#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"

#include <cstdint>
#include <vector>

using namespace std;
typedef ConcatenatedAlignedArrayRMM<float, float, std::uint64_t>::Quantity_t Quantity_t;

TEST_SUITE_REGISTRATION(ConcatenatedAlignedArrayRMMTest);

ConcatenatedAlignedArrayRMMTest::ConcatenatedAlignedArrayRMMTest() {
}

ConcatenatedAlignedArrayRMMTest::~ConcatenatedAlignedArrayRMMTest() {
}


void ConcatenatedAlignedArrayRMMTest::testAlignment() {

	Log::global_log->info() << "Testing AlignedArrayTriplet." << std::endl ;

	for (int i = 0; i < 11; ++i) {
		const size_t length = i;
		ConcatenatedAlignedArrayRMM<float, float, std::uint64_t> a(length);
		float * rx = a.begin_calc(Quantity_t::RX);
		float * ry = a.begin_calc(Quantity_t::RY);
		float * rz = a.begin_calc(Quantity_t::RZ);
		float * vx = a.begin_accum(Quantity_t::VX);
		float * vy = a.begin_accum(Quantity_t::VY);
		float * vz = a.begin_accum(Quantity_t::VZ);
		std::uint64_t * uid = a.begin_uid(Quantity_t::UID);

		ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(rx) % CACHE_LINE_SIZE), 0);
		ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(ry) % VCP_ALIGNMENT), 0);
		ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(rz) % VCP_ALIGNMENT), 0);
		ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(vx) % VCP_ALIGNMENT), 0);
		ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(vy) % VCP_ALIGNMENT), 0);
		ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(vz) % VCP_ALIGNMENT), 0);
		ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(uid) % VCP_ALIGNMENT), 0);
	}
}

void ConcatenatedAlignedArrayRMMTest::testGetters() {
	ConcatenatedAlignedArrayRMM<float, float, std::uint64_t> a(5);
	for (int i = static_cast<int>(Quantity_t::RX); i < static_cast<int>(Quantity_t::VX); ++i) {
		Quantity_t j = static_cast<Quantity_t>(i);
		a.get_calc(j,0) = i;
	}
	for (int i = static_cast<int>(Quantity_t::VX); i < static_cast<int>(Quantity_t::UID); ++i) {
		Quantity_t j = static_cast<Quantity_t>(i);
		a.get_accum(j,0) = i;
	}
	a.get_uid(Quantity_t::UID, 0) = 13;

	ASSERT_DOUBLES_EQUAL(0.0, a.get_calc(Quantity_t::RX, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(1.0, a.get_calc(Quantity_t::RY, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(2.0, a.get_calc(Quantity_t::RZ, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(3.0, a.get_accum(Quantity_t::VX, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(4.0, a.get_accum(Quantity_t::VY, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(5.0, a.get_accum(Quantity_t::VZ, 0), 1e-10);
	ASSERT_EQUAL(13ul, a.get_uid(Quantity_t::UID, 0));

	// try out the const methods
	const ConcatenatedAlignedArrayRMM<float, float, std::uint64_t> & b = a;
	ASSERT_DOUBLES_EQUAL(0.0, b.get_calc(Quantity_t::RX, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(1.0, b.get_calc(Quantity_t::RY, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(2.0, b.get_calc(Quantity_t::RZ, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(3.0, b.get_accum(Quantity_t::VX, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(4.0, b.get_accum(Quantity_t::VY, 0), 1e-10);
	ASSERT_DOUBLES_EQUAL(5.0, b.get_accum(Quantity_t::VZ, 0), 1e-10);
	ASSERT_EQUAL(13ul, b.get_uid(Quantity_t::UID, 0));
}

void ConcatenatedAlignedArrayRMMTest::testZero() {
	ConcatenatedAlignedArrayRMM<float, float, std::uint64_t> a(8);
	const int qbegin = static_cast<int>(Quantity_t::RX);
	const int qmid = static_cast<int>(Quantity_t::VX);
	const int qend = static_cast<int>(Quantity_t::UID);
	for (int q = qbegin; q < qmid; ++q) {
		Quantity_t qt1 = static_cast<Quantity_t>(q);
		for(size_t i = 0; i < 8; ++i) {
			a.get_calc(qt1, i)= 1.0;
		}
	}
	for (int q = qmid; q < qend; ++q) {
		Quantity_t qt1 = static_cast<Quantity_t>(q);
		for(size_t i = 0; i < 8; ++i) {
			a.get_accum(qt1, i)= 1.0;
		}
	}
	for(size_t i = 0; i < 8; ++i) {
		a.get_uid(Quantity_t::UID, i)= 1;
	}

	a.zero(5);
	for (int q = qbegin; q < qmid; ++q) {
		Quantity_t qt1 = static_cast<Quantity_t>(q);
		for(size_t i = 0; i < 5; ++i) {
			ASSERT_DOUBLES_EQUAL(1.0, a.get_calc(qt1, i), 1e-7);
		}
		for(size_t i = 5; i < 8; ++i) {
			ASSERT_DOUBLES_EQUAL(0.0, a.get_calc(qt1, i), 1e-7);
		}
	}
	for (int q = qmid; q < qend; ++q) {
		Quantity_t qt1 = static_cast<Quantity_t>(q);
		for(size_t i = 0; i < 5; ++i) {
			ASSERT_DOUBLES_EQUAL(1.0, a.get_accum(qt1, i), 1e-7);
		}
		for(size_t i = 5; i < 8; ++i) {
			ASSERT_DOUBLES_EQUAL(0.0, a.get_accum(qt1, i), 1e-7);
		}
	}
	for(size_t i = 0; i < 5; ++i) {
		ASSERT_EQUAL(1ul, a.get_uid(Quantity_t::UID, i));
	}
	for(size_t i = 5; i < 8; ++i) {
		ASSERT_EQUAL(0ul, a.get_uid(Quantity_t::UID, i));
	}

	a.zero(0);
	for (int q = qbegin; q < qmid; ++q) {
		Quantity_t qt1 = static_cast<Quantity_t>(q);
		for(size_t i = 0; i < 8; ++i) {
			ASSERT_DOUBLES_EQUAL(0.0, a.get_calc(qt1, i), 1e-7);
		}
	}
	for (int q = qmid; q < qend; ++q) {
		Quantity_t qt1 = static_cast<Quantity_t>(q);
		for(size_t i = 0; i < 8; ++i) {
			ASSERT_DOUBLES_EQUAL(0.0, a.get_accum(qt1, i), 1e-7);
		}
	}
	for(size_t i = 0; i < 8; ++i) {
		ASSERT_EQUAL(0ul, a.get_uid(Quantity_t::UID, i));
	}
}

template<typename TypeCalc, typename TypeAccum, typename TypeUID>
void check(const vector<TypeCalc>& x, const vector<TypeCalc>& y, const vector<TypeCalc>& z, const ConcatenatedAlignedArrayRMM<TypeCalc, TypeAccum, TypeUID>& A, int i) {

	for (int j = 0; j < i; ++j) {
		mardyn_assert(x[j] == A.get_calc(Quantity_t::RX,j));
		mardyn_assert(y[j] == A.get_calc(Quantity_t::RY,j));
		mardyn_assert(z[j] == A.get_calc(Quantity_t::RZ,j));

		mardyn_assert(x[j] == A.get_accum(Quantity_t::VX,j));
		mardyn_assert(y[j] == A.get_accum(Quantity_t::VY,j));
		mardyn_assert(z[j] == A.get_accum(Quantity_t::VZ,j));

		ASSERT_DOUBLES_EQUAL(x[j], A.get_calc(Quantity_t::RX,j), 0.0);
		ASSERT_DOUBLES_EQUAL(y[j], A.get_calc(Quantity_t::RY,j), 0.0);
		ASSERT_DOUBLES_EQUAL(z[j], A.get_calc(Quantity_t::RZ,j), 0.0);

		ASSERT_DOUBLES_EQUAL(x[j], A.get_accum(Quantity_t::VX,j), 0.0);
		ASSERT_DOUBLES_EQUAL(y[j], A.get_accum(Quantity_t::VY,j), 0.0);
		ASSERT_DOUBLES_EQUAL(z[j], A.get_accum(Quantity_t::VZ,j), 0.0);

	}

	for (int j = 0; j < i; ++j) {
		mardyn_assert(x[j] == A.get_uid(Quantity_t::UID,j));
	}
}

void ConcatenatedAlignedArrayRMMTest::testAppending() {
	vector<float> x, y, z;
	for(int i = 0; i < 17; ++i) {
		x.push_back(static_cast<float>(i));
		y.push_back(static_cast<float>(i+1));
		z.push_back(static_cast<float>(i+2));
	}

	ConcatenatedAlignedArrayRMM<float, float, std::uint64_t> A;
	for(int i = 0; i < 17; ++i) {
		std::array<float,3> val = {
			static_cast<float>(i),
			static_cast<float>(i+1),
			static_cast<float>(i+2)
		};
		std::array<float, 3> val2 = {
			static_cast<float>(i),
			static_cast<float>(i+1),
			static_cast<float>(i+2)
		};
		A.appendValues(val, val2, i, i);
	}

	float * rx = A.begin_calc(Quantity_t::RX);
	float * ry = A.begin_calc(Quantity_t::RY);
	float * rz = A.begin_calc(Quantity_t::RZ);
	float * vx = A.begin_accum(Quantity_t::VX);
	float * vy = A.begin_accum(Quantity_t::VY);
	float * vz = A.begin_accum(Quantity_t::VZ);
	std::uint64_t * uid = A.begin_uid(Quantity_t::UID);

	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(rx) % CACHE_LINE_SIZE), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(ry) % VCP_ALIGNMENT), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(rz) % VCP_ALIGNMENT), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(vx) % VCP_ALIGNMENT), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(vy) % VCP_ALIGNMENT), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(vz) % VCP_ALIGNMENT), 0);
	ASSERT_EQUAL(static_cast<int>(reinterpret_cast<intptr_t>(uid) % VCP_ALIGNMENT), 0);

	check(x,y,z,A,17);
}

void ConcatenatedAlignedArrayRMMTest::testIncreasingStorage() {
	vector<float> x, y, z;

	for(int i = 0; i < 111; ++i) {
		x.push_back(static_cast<float>(i));
		y.push_back(static_cast<float>(-i-1));
		z.push_back(static_cast<float>(i+2));
	}

	ConcatenatedAlignedArrayRMM<float, float, std::uint64_t> A;
	for(int i = 0; i < 11; ++i) {
		std::array<float, 3> val = {
			static_cast<float>(i),
			static_cast<float>(-i-1),
			static_cast<float>(i+2)
		};
		std::array<float, 3> val2 = {
			static_cast<float>(i),
			static_cast<float>(-i-1),
			static_cast<float>(i+2)
		};
		A.appendValues(val, val2, i, i);
		check(x, y, z, A, i);
	}

	A.increaseStorage(11, 2);

	for(int i = 11; i < 22; ++i) {
		std::array<float, 3> val = {
			static_cast<float>(i),
			static_cast<float>(-i-1),
			static_cast<float>(i+2)
		};
		std::array<float, 3> val2 = {
			static_cast<float>(i),
			static_cast<float>(-i-1),
			static_cast<float>(i+2)
		};
		A.appendValues(val, val2, i, i);
		check(x, y, z, A, i);
	}

	A.increaseStorage(22, 7);

	for(int i = 22; i < 111; ++i) {
		std::array<float, 3> val = {
			static_cast<float>(i),
			static_cast<float>(-i-1),
			static_cast<float>(i+2)
		};
		std::array<float, 3> val2 = {
			static_cast<float>(i),
			static_cast<float>(-i-1),
			static_cast<float>(i+2)
		};
		A.appendValues(val, val2, i, i);
		check(x, y, z, A, i);
	}

	check(x, y, z, A, 111);
}

