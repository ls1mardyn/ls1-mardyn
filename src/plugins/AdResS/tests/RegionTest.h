/*
 * Created on Fri Jul 26 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "utils/Testing.h"
#include "utils/TestWithSimulationSetup.h"

class RegionTest:public utils::TestWithSimulationSetup{
    TEST_SUITE(RegionTest);
    TEST_METHOD(testIsInnerPoint);
    TEST_METHOD(testIsBoxInHybrid);
    TEST_SUITE_END;

    public:
    RegionTest()=default;
    virtual ~RegionTest()=default;

    void testIsInnerPoint();

    
    void testIsBoxInHybrid();


};