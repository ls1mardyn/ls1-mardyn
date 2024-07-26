/*
 * Created on Fri Jul 26 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "utils/Testing.h"
#include "utils/TestWithSimulationSetup.h"
#include "plugins/AdResS/util/Region.h"

using namespace Resolution;

class RegionTest:public utils::Test{

    TEST_SUITE(RegionTest);
    TEST_METHOD(testIsInnerPoint);
    TEST_METHOD(testIsBoxInHybrid);
    TEST_METHOD(testComputeIntersection);
    TEST_SUITE_END;

    public:
    RegionTest();
    virtual ~RegionTest()=default;

    void testIsInnerPoint();
    void testIsBoxInHybrid();
    void testComputeIntersection();
    //void testIsRegionInBox();

    private:

    FPRegion region;
};