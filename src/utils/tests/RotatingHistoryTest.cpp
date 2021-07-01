#include "RotatingHistoryTest.h"
#include "../RotatingHistory.h"

#include <initializer_list>
#include <iostream>

TEST_SUITE_REGISTRATION(RotatingHistoryTest);

template <typename RH>
void testExpectation(RH&& rh, std::initializer_list<double>&& values) {
	size_t count{0ul};
	for (auto&& val : rh) {
		ASSERT_EQUAL(std::data(values)[count], val);
		++count;
	}
	size_t expected_count = values.size();
	ASSERT_EQUAL(expected_count, count);
}

void RotatingHistoryTest::testDefaultConstructed() {
	RotatingHistory<double> rh;

	rh.insert(1.);
	testExpectation(rh, {1.});
	testExpectation(std::as_const(rh), {1.});
}

void RotatingHistoryTest::testZeroSize() {
	RotatingHistory<double> rh(0ul);

	rh.insert(1.);
	testExpectation(rh, {});
	testExpectation(std::as_const(rh), {});
}

void RotatingHistoryTest::testFiveElements() {
	RotatingHistory<double> rh(5);

	rh.insert(1.);
	testExpectation(rh, {1.});
	testExpectation(std::as_const(rh), {1.});

	rh.insert(2.);
	rh.insert(4.);
	testExpectation(rh, {1., 2., 4.});
	testExpectation(std::as_const(rh), {1., 2., 4.});

	rh.insert(8.);
	rh.insert(16.);
	testExpectation(rh, {1., 2., 4., 8., 16.});
	testExpectation(std::as_const(rh), {1., 2., 4., 8., 16.});

	rh.insert(32.);
	rh.insert(64.);
	testExpectation(rh, {4., 8., 16., 32., 64.});
	testExpectation(std::as_const(rh), {4., 8., 16., 32., 64.});
}
