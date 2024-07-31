//
// Created by alex on 7/26/24.
//

#include "AdResSAveragerTest.h"

TEST_SUITE_REGISTRATION(AdResSAveragerTest);

std::size_t AdResSAveragerTest::MAX_STEPS = 10;
std::size_t AdResSAveragerTest::DATA_SIZE = 16;

AdResSAveragerTest::AdResSAveragerTest() {
    std::vector<double> tmp (DATA_SIZE, 0.0);
    _avg.setDataSize(tmp);

    for (int i = 0; i < DATA_SIZE; i++) tmp[i] = i + 1;
    for (int j = 0; j < MAX_STEPS; j++) _avg.averageData(tmp);
}

void AdResSAveragerTest::checkCount() {
    ASSERT_EQUAL((int) MAX_STEPS, _avg.getStepCount());
}

void AdResSAveragerTest::checkSize() {
    ASSERT_EQUAL((int) DATA_SIZE, static_cast<int>(_avg.getSumData().size()));
}

void AdResSAveragerTest::checkAverage() {
    std::vector<double> data (DATA_SIZE, 0.0);
    _avg.getAveragedData(data);

    for (int idx = 0; idx < DATA_SIZE; idx++) {
        double target = idx + 1;
        ASSERT_DOUBLES_EQUAL(target, data[idx], 1e-15);
    }
}

void AdResSAveragerTest::checkAverageCopy() {
    auto data = _avg.getAveragedDataCopy();

    for (int idx = 0; idx < DATA_SIZE; idx++) {
        double target = idx + 1;
        ASSERT_DOUBLES_EQUAL(target, data[idx], 1e-15);
    }
}

void AdResSAveragerTest::checkSum() {
    const auto& sumData = _avg.getSumData();
    for (int idx = 0; idx < DATA_SIZE; idx++) {
        double target = MAX_STEPS * (idx + 1);
        ASSERT_EQUAL(target, sumData[idx]);
    }
}

void AdResSAveragerTest::checkReset() {
    _avg.reset();
    const auto& sumData = _avg.getSumData();
    for (int idx = 0; idx < DATA_SIZE; idx++) {
        ASSERT_EQUAL(0.0, sumData[idx]);
    }
    ASSERT_EQUAL(0, _avg.getStepCount());
}
