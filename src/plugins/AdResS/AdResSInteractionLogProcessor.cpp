//
// Created by alex on 09.05.23.
//

#include "AdResSInteractionLogProcessor.h"

AdResSInteractionLogProcessor::AdResSInteractionLogProcessor(const double cutoffRadius, const double ljCutoffRadius,
                                                             const vector<FPRegion> &fpRegions) : CellProcessor(
        cutoffRadius, ljCutoffRadius), _fpRegions(fpRegions), _processed(false) { }
AdResSInteractionLogProcessor::~AdResSInteractionLogProcessor() = default;
void AdResSInteractionLogProcessor::initTraversal() {}
void AdResSInteractionLogProcessor::preprocessCell(ParticleCell &cell) {}
double AdResSInteractionLogProcessor::processSingleMolecule(Molecule *m1, ParticleCell &cell2) {return 0;}
void AdResSInteractionLogProcessor::postprocessCell(ParticleCell &cell) {}

void AdResSInteractionLogProcessor::endTraversal() { _processed = true; }

static inline bool overlapped(std::array<double, 3> low1, std::array<double, 3> high1, std::array<double, 3> low2, std::array<double, 3> high2) {
    return low1[0] <= high2[0] && low2[0] <= high1[0] &&
           low1[1] <= high2[1] && low2[1] <= high1[1] &&
           low1[2] <= high2[2] && low2[2] <= high1[2];
}

void AdResSInteractionLogProcessor::processCell(ParticleCell &cell) {
    auto low = cell.getBoxMinArray();
    auto high = cell.getBoxMaxArray();

    for(unsigned long index = 0; index < _fpRegions.size(); index++) {
        const auto& fpr = _fpRegions[index];
        //check for overlapping boxes
        if(overlapped(fpr._lowHybrid, fpr._highHybrid, low, high))
            _fprID_to_cells[index].emplace(&cell);
    }
}


void AdResSInteractionLogProcessor::processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) {
    auto low1 = cell1.getBoxMinArray();
    auto high1 = cell1.getBoxMaxArray();
    auto low2 = cell2.getBoxMinArray();
    auto high2 = cell2.getBoxMaxArray();

    for(unsigned long index = 0; index < _fpRegions.size(); index++) {
        const auto& fpr = _fpRegions[index];
        //check for overlapping boxes
        bool ov1, ov2;
        ov1 = overlapped(fpr._lowHybrid, fpr._highHybrid, low1, high1);
        ov2 = overlapped(fpr._lowHybrid, fpr._highHybrid, low2, high2);
        if(ov1) {
            _fprID_to_cells[index].emplace(&cell1);
        }
        if(ov2) {
            _fprID_to_cells[index].emplace(&cell2);
        }
        if(ov1 || ov2) {
            _fprID_to_cellpairs[index].emplace(&cell1, &cell2);
        }

    }
}
