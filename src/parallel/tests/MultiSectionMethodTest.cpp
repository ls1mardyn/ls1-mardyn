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

void MultiSectionMethodTest::testGetCoordsFromRank() {
	int numProcs = 2 * 3 * 5;
	std::array<size_t, 3> grid{2, 3, 5};
	for (int rank = 0; rank < numProcs; rank++) {
		auto coord = MultiSectionMethod::getCoordsFromRank(grid, rank);
		auto shouldBeRank = coord[0] * grid[1] * grid[2] + coord[1] * grid[2] + coord[2];
		ASSERT_TRUE_MSG("Coords should match", shouldBeRank == rank);
	}
}
void MultiSectionMethodTest::testInitializeRegularGrid() {
	std::array<double, 3> domainLength{.1, 3.5, .7}, boxMin{}, boxMax{};
	std::array<size_t, 3> gridSize{3, 7, 9}, gridCoords{2, 4, 6};
	std::tie(boxMin, boxMax) = MultiSectionMethod::initializeRegularGrid(domainLength, gridSize, gridCoords);
	ASSERT_DOUBLES_EQUAL(domainLength[0] * 2. / 3., boxMin[0], 1e-15);
	ASSERT_EQUAL(domainLength[0], boxMax[0]);
	ASSERT_DOUBLES_EQUAL(domainLength[1] * 4. / 7., boxMin[1], 1e-15);
	ASSERT_DOUBLES_EQUAL(domainLength[1] * 5. / 7., boxMax[1], 1e-15);
	ASSERT_DOUBLES_EQUAL(domainLength[2] * 6. / 9., boxMin[2], 1e-15);
	ASSERT_DOUBLES_EQUAL(domainLength[2] * 7. / 9., boxMax[2], 1e-15);
}
