//
// Created by alex on 7/26/24.
//

#include "AdResSWeightFunctionTest.h"

TEST_SUITE_REGISTRATION(AdResSWeightFunctionTest);

AdResSWeightFunctionTest::AdResSWeightFunctionTest() {
    _region = Resolution::FPRegion({2, 2, 2}, {3, 3, 3}, {1, 1, 1});
    _region.init();
}

void AdResSWeightFunctionTest::testBounds() {
    for (auto f : Weight::functions) {
        if (f == Weight::flat) break;
        checkBounds(f);
    }
}

void AdResSWeightFunctionTest::testCenter() {
    for (auto f : Weight::functions) {
        if (f == Weight::flat) break;
        checkCenter(f);
    }
}

void AdResSWeightFunctionTest::testOther() {
    for (auto f : Weight::functions) {
        if (f == Weight::flat) break;
        checkOther(f);
    }
}

void AdResSWeightFunctionTest::checkBounds(Weight::function_t fun) {
    // FP inter C inter CG
    for (const auto& point_set : _test_sets) {
        double wFP = fun(point_set[0], _region);
        double wCG = fun(point_set[4], _region);
        ASSERT_EQUAL(wFP, 1.0);
        ASSERT_EQUAL(wCG, 0.0);
    }
}

void AdResSWeightFunctionTest::checkCenter(Weight::function_t fun) {
    // FP inter C inter CG
    for (const auto& point_set : _test_sets) {
        double wH = fun(point_set[2], _region);
        ASSERT_DOUBLES_EQUAL(0.5, wH, 0.3);
    }
}

void AdResSWeightFunctionTest::checkOther(Weight::function_t fun) {
    // FP inter C inter CG
    for (const auto& point_set : _test_sets) {
        double wFP = fun(point_set[0], _region);
        double wFPH = fun(point_set[1], _region);
        double wH = fun(point_set[2], _region);
        double wHCG = fun(point_set[3], _region);
        double wCG = fun(point_set[4], _region);
        ASSERT_TRUE(wFP >= wFPH);
        ASSERT_TRUE(wFPH >= wH);
        ASSERT_TRUE(wH >= wHCG);
        ASSERT_TRUE(wHCG >= wCG);
    }
}
