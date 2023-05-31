//
// Created by alex on 31.05.23.
//

#ifndef MARDYN_ADRESSTEST_H
#define MARDYN_ADRESSTEST_H

#include "utils/TestWithSimulationSetup.h"

class AdResSTest : public utils::TestWithSimulationSetup {
TEST_SUITE(AdResSTest);
        TEST_METHOD(computeForcesTest);
    TEST_SUITE_END;

public:
    AdResSTest();
    virtual ~AdResSTest();

    /**
     * Checks if all forces are computed correctly according to AdResS scheme.
     * */
    void computeForcesTest();
};


#endif //MARDYN_ADRESSTEST_H
