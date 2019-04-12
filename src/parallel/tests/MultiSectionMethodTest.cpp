/**
 * @file MultiSectionMethodTest.cpp
 * @author seckler
 * @date 12.04.19
 */

#include "MultiSectionMethodTest.h"
#include "parallel/MultiSectionMethod.h"

TEST_SUITE_REGISTRATION(MultiSectionMethodTest);

void MultiSectionMethodTest::testGetOptimalGrid() {
	int numProcs = 2 * 3 * 5;
	testGetOptimalGridBody({1., 3., 5.}, numProcs);
	testGetOptimalGridBody({1., 5., 3.}, numProcs);
	testGetOptimalGridBody({3., 1., 5.}, numProcs);
	testGetOptimalGridBody({3., 5., 1.}, numProcs);
	testGetOptimalGridBody({5., 1., 3.}, numProcs);
	testGetOptimalGridBody({5., 3., 1.}, numProcs);
}
void MultiSectionMethodTest::testGetOptimalGridBody(const std::array<double, 3> domainLength, int numProcs) {
	auto grid = MultiSectionMethod::getOptimalGrid(domainLength, numProcs);

	ASSERT_TRUE_MSG("ordering wrong", (domainLength[0] < domainLength[1]) == (grid[0] < grid[1]));
	ASSERT_TRUE_MSG("ordering wrong", (domainLength[0] < domainLength[2]) == (grid[0] < grid[2]));
	ASSERT_TRUE_MSG("ordering wrong", (domainLength[1] < domainLength[2]) == (grid[1] < grid[2]));
}
