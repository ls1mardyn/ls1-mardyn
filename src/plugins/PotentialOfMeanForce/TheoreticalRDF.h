//
// Created by alex on 3/12/25.
//

#ifndef MARDYN_THEORETICALRDF_H
#define MARDYN_THEORETICALRDF_H

#include <vector>
#include <array>

#include "particleContainer/ParticleContainer.h"

class TheoreticalRDF {
public:
    TheoreticalRDF();

    TheoreticalRDF(const std::array<double, 3>& low, const std::array<double, 3>& high, int bins);

    void measure(ParticleContainer* pc);

    std::vector<double> getCurrentRDF();

    void writeFile(const std::string& name, unsigned long simstep);

private:
    std::array<double, 3> regionLow;
    std::array<double, 3> regionHigh;
    std::array<double, 3> outerLow;
    std::array<double, 3> outerHigh;
    int numBins;
    double binSize;
    double rhoGlobal;

    std::vector<int> counts;
    std::vector<std::vector<int>> threadCounts;
    unsigned long countedParticles;
    double cutoff2;
};


#endif //MARDYN_THEORETICALRDF_H
