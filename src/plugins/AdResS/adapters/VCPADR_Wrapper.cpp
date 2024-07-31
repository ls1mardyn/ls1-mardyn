//
// Created by alex on 7/15/24.
//

#include "VCPADR_Wrapper.h"
#include "../features/Resolution.h"
#include "particleContainer/FullParticleCell.h"
#include "Simulation.h"

VCPADR_Wrapper::VCPADR_Wrapper(Domain &domain, double cutoffRadius, double LJcutoffRadius,
                               const Resolution::Handler &resolutionHandler) :
    CellProcessor(cutoffRadius, LJcutoffRadius),
    _resolution_handler(resolutionHandler),
    _is_init(false) {}

void VCPADR_Wrapper::init() {
    _reference_processor = std::make_unique<VectorizedCellProcessor>(*_simulation.getDomain(), getCutoffRadius(), getLJCutoffRadius());
    _adr_processor = std::make_unique<VCPADR>(*_simulation.getDomain(), getCutoffRadius(), getLJCutoffRadius(), _resolution_handler);
    dynamic_cast<VCPADR*>(_adr_processor.get())->init();

    int numThreads = mardyn_get_max_threads();;
    _cell_maps.resize(numThreads);
}

void VCPADR_Wrapper::initTraversal() {
    _adr_processor->initTraversal();
    _reference_processor->initTraversal();
}

void VCPADR_Wrapper::processCell(ParticleCell &cell) {
    if (checkCell(cell)) _adr_processor->processCell(cell);
    else _reference_processor->processCell(cell);
}

void VCPADR_Wrapper::processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) {
    if (checkCell(cell1) || checkCell(cell2)) _adr_processor->processCellPair(cell1, cell2, sumAll);
    else _reference_processor->processCellPair(cell1, cell2, sumAll);
}

void VCPADR_Wrapper::endTraversal() {
    _adr_processor->endTraversal();
    _reference_processor->endTraversal();

    if (_is_init) return;
    // at this point all cells should have been checked by at least one thread
    // merge all thread local maps, to avoid issues when using dynamic scheduling or blocking

    // first merge all into first buffer
    for (int tid = 1; tid < mardyn_get_max_threads(); tid++) {
        for (auto [cell, value] : _cell_maps[tid]) {
            _cell_maps[0][cell] = value;
        }
    }
    // broadcast
    for (int tid = 1; tid < mardyn_get_max_threads(); tid++) {
        _cell_maps[tid] = _cell_maps[0];
    }

    _is_init = true;
}

bool VCPADR_Wrapper::checkCell(ParticleCell &cell) {
    int tid = mardyn_get_thread_num();

    if (_is_init) return _cell_maps[tid][&cell];

    // unknown -> insert
    if (_cell_maps[tid].find(&cell) == _cell_maps[tid].end()) {
        auto& region = _resolution_handler.getRegions()[0];
        _cell_maps[tid][&cell] = region.isBoxInHybrid(cell.getBoxMinArray(), cell.getBoxMaxArray());
    }

    return _cell_maps[tid][&cell];
}

