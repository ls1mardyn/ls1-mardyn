#include "RegionTest.h"

TEST_SUITE_REGISTRATION(RegionTest);

void RegionTest::testIsInnerPoint(){
    std::array<double, 3> high = {10,10,10};
    std::array<double, 3> low = {0,0,0};
    std::array<double, 3> dims = {5,5,5}; 
    std::array<double, 3> pos = {1,1,1};
    FPRegion region{low,high,dims};
    region.init();

    region.isInnerPoint(pos,low,high);

    ASSERT_TRUE(region.isInnerPoint(pos,low,high));

}

void RegionTest::testIsInsideResolutionRegion(){

}

void RegionTest::testParticleInHalo(){

}