/*
 * ExtendedProfileSamplingTest.cpp
 *
 *  Created on: Jun 2022
 *      Author: homes
 */

#ifndef DEXTENDEDPROFILESAMPLINGTEST_H
#define DEXTENDEDPROFILESAMPLINGTEST_H

#include "utils/TestWithSimulationSetup.h"
#include "plugins/ExtendedProfileSampling.h"

// Test of ExtendedProfileSampling plugin
// correct sampling of only a few quantities is tested as other quantities may require special setting, e.g. legacy cell processor

class ExtendedProfileSamplingTest : public utils::TestWithSimulationSetup {

    TEST_SUITE(ExtendedProfileSamplingTest);
    TEST_METHOD(testEPSampling);
    TEST_SUITE_END;

public:

    ExtendedProfileSamplingTest();

    virtual ~ExtendedProfileSamplingTest();

    void testEPSampling();

};

#endif //DEXTENDEDPROFILESAMPLINGTEST_H
