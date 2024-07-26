#include "RegionTest.h"
#include "plugins/AdResS/util/Region.h"
using namespace Resolution;

TEST_SUITE_REGISTRATION(RegionTest);

void RegionTest::testIsInnerPoint(){

    std::array<double, 3> low = {0,0,0};    
    std::array<double, 3> high = {10,10,10};
    std::array<double, 3> dims = {0,0,0};
    std::array<double, 3> p1 = {5,5,5};
    std::array<double, 3> p2 = {11,11,11};

    FPRegion region{low,high,dims};
    region.init();

    ASSERT_TRUE(region.isInnerPoint(p1,low,high));
    ASSERT_TRUE(!region.isInnerPoint(p2,low,high));

}

void RegionTest::testIsBoxInHybrid(){
    std::array<double, 3> low = {0,0,0};    
    std::array<double, 3> high = {10,10,10};
    std::array<double, 3> dims = {5,5,5};
    FPRegion region{low,high,dims};
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