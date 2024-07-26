#include "RegionTest.h"


TEST_SUITE_REGISTRATION(RegionTest);

RegionTest::RegionTest(){
    region._low = {0,0,0};
    region._high = {10,10,10};
    region._hybridDims = {0,0,0};
    region.init();
}

void RegionTest::testComputeIntersection(){
    //TODO: think if tests are needed.
}

void RegionTest::testIsInnerPoint(){
    
    std::array<double, 3> p1 = {5,5,5};
    std::array<double, 3> p2 = {11,11,11};


    ASSERT_TRUE(region.isInnerPoint(p1,region._low,region._high));
    ASSERT_TRUE(!region.isInnerPoint(p2,region._low,region._high));
}

void RegionTest::testIsBoxInHybrid(){
    region._hybridDims = {5,5,5};
    region.init();


    std::array<double, 3> cell_low_out = {-20,-20,-20};
    std::array<double, 3> cell_high_out = {-15,-15,-15};

    ASSERT_TRUE(!region.isBoxInHybrid(cell_low_out,cell_high_out));

    std::array<double, 3> cell_low_in = {10,10,10};
    std::array<double, 3> cell_high_in = {12,12,-12};

    ASSERT_TRUE(region.isBoxInHybrid(cell_low_in,cell_high_in));

    std::array<double, 3> cell_low_partial = {-6,-6,-6};
    std::array<double, 3> cell_high_partial = {-3,-3,-3};

    ASSERT_TRUE(region.isBoxInHybrid(cell_low_partial,cell_high_partial));

    std::array<double, 3> cell_low_larger = {-20,-20,-20};
    std::array<double, 3> cell_high_larger = {20,20,20};

    ASSERT_TRUE_MSG("cell cannot contain hybrid resolution region",!region.isBoxInHybrid(cell_low_larger,cell_high_larger));



}