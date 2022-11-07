#pragma once

#include "utils/TestWithSimulationSetup.h"
#include "plugins/NEMD/VelocityExchange.h"

// Test of VelocityExchange plugin
// Exchange of velocities of coldest particle in warmest region with warmest particle in coldest region is tested

class VelocityExchangeTest : public utils::TestWithSimulationSetup {

    TEST_SUITE(VelocityExchangeTest);
    TEST_METHOD(testExchangeVelocities);
    TEST_SUITE_END;

public:

    VelocityExchangeTest();

    virtual ~VelocityExchangeTest();

    void testExchangeVelocities();

};
