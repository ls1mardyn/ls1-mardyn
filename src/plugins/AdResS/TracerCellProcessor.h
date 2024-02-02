//
// Created by alex on 31.01.24.
//

#ifndef MARDYN_TRACERCELLPROCESSOR_H
#define MARDYN_TRACERCELLPROCESSOR_H

#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "AdResSData.h"

class TracerCellProcessor : public LegacyCellProcessor {
public:
    TracerCellProcessor(const double cutoffRadius, const double ljCutoffRadius,
                        ParticlePairsHandler *particlePairsHandler, std::vector<Resolution>& ctr);

private:
    void processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) override;

    double processSingleMolecule(Molecule *m1, ParticleCell &cell2) override;

    void processCell(ParticleCell &cell) override;

    std::vector<Resolution>& _comp_to_res;
    ParticlePairsHandler* _particlePairsHandler;
};


#endif //MARDYN_TRACERCELLPROCESSOR_H
