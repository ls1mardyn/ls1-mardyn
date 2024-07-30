/*
 * Created on Tue Jul 30 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "utils/Testing.h"
#include "utils/TestWithSimulationSetup.h"
#include "plugins/AdResS/util/Region.h"
#include "plugins/AdResS/features/Resolution.h"
#include "plugins/AdResS/AdResS.h"

#include "particleContainer/TraversalTuner.h"

using namespace Resolution;

class ResolutionTest: public utils::TestWithSimulationSetup{
    TEST_SUITE(ResolutionTest);
    TEST_METHOD(testCheckResolution);
    TEST_SUITE_END;

    public:
    ResolutionTest();
    virtual ~ResolutionTest()=default;

    void testCheckResolution();

    private:
    std::string file_name = "AdResS-empty-10x10x10.inp";

    Handler handler;
    Config config;
    
};