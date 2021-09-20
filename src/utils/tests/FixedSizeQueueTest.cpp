#include "FixedSizeQueueTest.h"
#include "../FixedSizeQueue.h"

#include <initializer_list>
#include <iostream>

TEST_SUITE_REGISTRATION(FixedSizeQueueTest);

template <typename RH>
void testExpectationHelper(RH&& rh, const std::initializer_list<double>& values) {
	size_t count{0ul};
	for (auto&& val : rh) {
		ASSERT_EQUAL(std::data(values)[count], val);
		++count;
	}
	size_t expected_count = values.size();
	ASSERT_EQUAL(expected_count, count);
}

template <typename RH>
void testExpectation(RH&& rh, std::initializer_list<double>&& values) {
	testExpectationHelper(rh, values);
	testExpectationHelper(std::as_const(rh), values);
}

void FixedSizeQueueTest::testDefaultConstructed() {
	FixedSizeQueue<double> rh;
	// Should have capacity 1!

	testExpectation(rh, {});

	rh.insert(1.);
	testExpectation(rh, {1.});

	rh.insert(2.);
	testExpectation(rh, {2.});
}

void FixedSizeQueueTest::testZeroSize() {
	FixedSizeQueue<double> rh(0ul);

	rh.insert(1.);
	testExpectation(rh, {});
}

void FixedSizeQueueTest::testFiveElements() {
	FixedSizeQueue<double> rh(5);

	rh.insert(1.);
	testExpectation(rh, {1.});

	rh.insert(2.);
	rh.insert(4.);
	testExpectation(rh, {1., 2., 4.});

	rh.insert(8.);
	rh.insert(16.);
	testExpectation(rh, {1., 2., 4., 8., 16.});

	rh.insert(32.);
	rh.insert(64.);
	testExpectation(rh, {4., 8., 16., 32., 64.});
}

void FixedSizeQueueTest::testGrowth() {
	FixedSizeQueue<double> rh(1);

	rh.insert(1.);
	testExpectation(rh, {1.});

	rh.insert(2.);
	testExpectation(rh, {2.});

	rh.setCapacity(3);
	rh.insert(16.);
	rh.insert(32.);
	testExpectation(rh, {2., 16., 32.});

	rh.insert(64.);
	testExpectation(rh, {16., 32., 64.});
}

void FixedSizeQueueTest::testShrink() {
	FixedSizeQueue<double> rh(4);

	rh.insert(1.);
	rh.insert(2.);
	rh.insert(3.);
	rh.insert(4.);
	testExpectation(rh, {1., 2., 3., 4.});

	rh.setCapacity(1ul);
	testExpectation(rh, {4.});

	rh.insert(16.);
	testExpectation(rh, {16.});
}