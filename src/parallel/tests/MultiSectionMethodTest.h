/**
 * @file MultiSectionMethodTest.h
 * @author seckler
 * @date 12.04.19
 */

#pragma once
#include "utils/TestWithSimulationSetup.h"

class MultiSectionMethodTest : public utils::TestWithSimulationSetup {
	TEST_SUITE(MultiSectionMethodTest);
	TEST_METHOD(testGetOptimalGrid);
	TEST_METHOD(testGetCoordsFromRank);
	TEST_METHOD(testInitializeRegularGrid);
	TEST_SUITE_END();

public:
	void testGetOptimalGrid();
	void testGetCoordsFromRank();
	void testInitializeRegularGrid();

private:
	void testGetOptimalGridBody(std::array<double, 3> domainLength, int numProcs);
};
