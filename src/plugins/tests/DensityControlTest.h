#pragma once

#include "utils/TestWithSimulationSetup.h"
#include "utils/Testing.h"
#include "plugins/NEMD/DensityControl.h"

// Test of DensityControl plugin
// Set of density within test region is tested

class DensityControlTest : public utils::TestWithSimulationSetup {

    TEST_SUITE(DensityControlTest);
    TEST_METHOD(testDensityControl);
    TEST_SUITE_END;

public:

    DensityControlTest();

    virtual ~DensityControlTest();

    void testDensityControl();

};
