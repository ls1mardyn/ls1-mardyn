//
// Created by kruegener on 5/31/2018.
//


#ifndef DCOMALIGNERTEST_H
#define DCOMALIGNERTEST_H

#include "utils/TestWithSimulationSetup.h"
#include "plugins/COMaligner.h"

class COMalignerTest : public utils::TestWithSimulationSetup {

    TEST_SUITE(COMalignerTest);
    TEST_METHOD(testCOMalign);
    TEST_SUITE_END;

public:

    COMalignerTest();

    virtual ~COMalignerTest();

    void testCOMalign();

};

#endif //DCOMALIGNERTEST_H
