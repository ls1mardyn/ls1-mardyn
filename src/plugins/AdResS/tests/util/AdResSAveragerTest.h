//
// Created by alex on 7/26/24.
//

#ifndef MARDYN_ADRESSAVERAGERTEST_H
#define MARDYN_ADRESSAVERAGERTEST_H

#include "utils/TestWithSimulationSetup.h"
#include "plugins/AdResS/util/Averager.h"

#include <vector>

class AdResSAveragerTest : public utils::Test {
TEST_SUITE(AdResSAveragerTest);
    TEST_METHOD(checkSize);
    TEST_METHOD(checkCount);
    TEST_METHOD(checkAverage);
    TEST_METHOD(checkAverageCopy);
    TEST_METHOD(checkSum);
    TEST_METHOD(checkReset);
TEST_SUITE_END;

public:
    AdResSAveragerTest();
    virtual ~AdResSAveragerTest() = default;

    void checkSize();
    void checkCount();
    void checkAverage();
    void checkAverageCopy();
    void checkSum();
    void checkReset();

private:
    Averager<std::vector<double>> _avg;
    static std::size_t DATA_SIZE;
    static std::size_t MAX_STEPS;
};

#endif //MARDYN_ADRESSAVERAGERTEST_H
