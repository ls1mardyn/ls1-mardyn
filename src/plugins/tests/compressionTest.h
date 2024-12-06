//
// Created by fernanor on 2018-12-13.
//


#ifndef COMPRESSIONTEST_H
#define COMPRESSIONTEST_H

#include <climits>
#include <exception>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>

#include "utils/TestWithSimulationSetup.h"
#include "plugins/compression.h"

class compressionTest : public utils::TestWithSimulationSetup {

    TEST_SUITE(compressionTest);
    TEST_METHOD(testNone);
#ifdef ENABLE_LZ4
    TEST_METHOD(testLz4Random);
    TEST_METHOD(testLz4Sine);
#endif
    TEST_METHOD(testThrowOnFailedCreate);
    TEST_SUITE_END;

public:

    compressionTest();

    virtual ~compressionTest();
    void testNone();
#ifdef ENABLE_LZ4
    void testLz4Random();
    void testLz4Sine();
#endif
    void testThrowOnFailedCreate();
};

#endif //COMPRESSIONTEST_H
