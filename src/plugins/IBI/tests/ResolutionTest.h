/*
 * Created on Thu Jul 25 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "utils/Testing.h"
#include "utils/TestWithSimulationSetup.h"

#include "../Resolution.h"


class ResolutionTest:public utils::TestWithSimulationSetup{
    TEST_SUITE(ResolutionTest);
    TEST_METHOD(testCheckResolution);
    TEST_SUITE_END;

    public:

    ResolutionTest()=default;
    virtual ~ResolutionTest()=default;

    void testCheckResolution();

};