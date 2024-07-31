//
// Created by alex on 7/30/24.
//

#ifndef MARDYN_VCPADR_WRAPPERTEST_H
#define MARDYN_VCPADR_WRAPPERTEST_H

#include "utils/TestWithSimulationSetup.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "WrapOpenMP.h"

class CellProcMock : public CellProcessor {
public:
    CellProcMock(const double cutoffRadius, const double ljCutoffRadius) : CellProcessor(cutoffRadius, ljCutoffRadius) {};
    ~CellProcMock() override {};
    void initTraversal() override {
        _thread_cells.clear();
        _thread_cellPairs.clear();
        _cells.clear();
        _cellPairs.clear();

        int numThreads = mardyn_get_max_threads();
        _thread_cells.resize(numThreads);
        _thread_cellPairs.resize(numThreads);
    };
    void endTraversal() override {
        for (auto& vec : _thread_cells) {
            for (auto* ptr : vec) {
                _cells.push_back(ptr);
            }
        }
        for (auto& vec : _thread_cellPairs) {
            for (auto& arr : vec) {
                _cellPairs.push_back(arr);
            }
        }
    };

    void preprocessCell(ParticleCell &cell) override {};
    double processSingleMolecule(Molecule *m1, ParticleCell &cell2) override { return 0; };
    void postprocessCell(ParticleCell &cell) override {};

    void processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) override {
        int tid = mardyn_get_thread_num();
        _thread_cellPairs[tid].push_back({&cell1, &cell2});
    };
    void processCell(ParticleCell &cell) override {
        int tid = mardyn_get_thread_num();
        _thread_cells[tid].push_back(&cell);
    };

    std::vector<ParticleCell*>& get_cells() { return _cells; }
    std::vector<std::array<ParticleCell*, 2>>& get_cellPairs() { return _cellPairs; }

private:
    std::vector<ParticleCell*> _cells;
    std::vector<std::vector<ParticleCell*>> _thread_cells;
    std::vector<std::array<ParticleCell*, 2>> _cellPairs;
    std::vector<std::vector<std::array<ParticleCell*, 2>>> _thread_cellPairs;
};

class VCPADR_WrapperTest : public utils::TestWithSimulationSetup {
// We only want to check if the correct CellProcessors are being used
// Force computation is checked for each CellProcessor individually
TEST_SUITE(VCPADR_WrapperTest);
    TEST_METHOD(processorSelectionTest);
TEST_SUITE_END;

public:
    VCPADR_WrapperTest() = default;
    virtual ~VCPADR_WrapperTest() = default;

    void processorSelectionTest();
};


#endif //MARDYN_VCPADR_WRAPPERTEST_H
