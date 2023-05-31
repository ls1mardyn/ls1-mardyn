//
// Created by alex on 31.05.23.
//

#ifndef MARDYN_ADRESSFORCEADAPTERTEST_H
#define MARDYN_ADRESSFORCEADAPTERTEST_H


#include "utils/TestWithSimulationSetup.h"

/**
 * Tests for the AdResSForceAdapter
 * Checks if molecule pairs are handled correctly regarding AdResS forces.
 * */
class AdResSForceAdapterTest : public utils::TestWithSimulationSetup {
TEST_SUITE(AdResSForceAdapterTest);
        TEST_METHOD(processPairTest);
    TEST_SUITE_END;

public:
    AdResSForceAdapterTest();
    virtual ~AdResSForceAdapterTest();

    /**
     * Creates a domain of 10x10x10 with rCutoff = 2.
     * Places molecules in all region (FP,H,CG) combinations and lets them interact with each other.
     * */
    void processPairTest();
};


#endif //MARDYN_ADRESSFORCEADAPTERTEST_H
