//
// Created by alex on 3/12/25.
//

#include "TheoreticalRDF.h"
#include "Simulation.h"
#include "Domain.h"

TheoreticalRDF::TheoreticalRDF() : TheoreticalRDF({0, 0, 0}, {0, 0, 0}, 0) { }

TheoreticalRDF::TheoreticalRDF(const std::array<double, 3> &low, const std::array<double, 3> &high, int bins) :
        regionLow(low), regionHigh(high), outerLow(low), outerHigh(high),
        numBins(bins), binSize(0), rhoGlobal(0), counts(bins, 0), threadCounts(mardyn_get_max_threads()), countedParticles(0), cutoff2(0) {
    rhoGlobal = _simulation.getDomain()->getglobalRho();
    double cutoff = _simulation.getcutoffRadius();
    cutoff2 = std::pow(cutoff, 2);
    binSize = cutoff / numBins;

    for (int i = 0; i < 3; i++) {
        if (regionLow[i] < 0 || regionHigh[i] > _simulation.getDomain()->getGlobalLength(i)) throw std::runtime_error("Region is out of bounds");

        outerLow[i] -= cutoff;
        outerHigh[i] += cutoff;
    }

    for (auto& buffer : threadCounts) buffer.resize(numBins, 0);
}

void TheoreticalRDF::measure(ParticleContainer *pc) {
    for (auto& buffer : threadCounts) std::fill(buffer.begin(), buffer.end(), 0);


    for (auto itRegion = pc->regionIterator(regionLow.data(), regionHigh.data(), ParticleIterator::ONLY_INNER_AND_BOUNDARY); itRegion.isValid(); ++itRegion) {
        #pragma omp parallel
        for (auto itOuter = pc->regionIterator(outerLow.data(), outerHigh.data(), ParticleIterator::ALL_CELLS); itOuter.isValid(); ++itOuter) {
            if (itRegion->getID() == itOuter->getID()) continue;
            auto& molRegion = *itRegion;
            auto& molOuter = *itOuter;

            const auto rRegion = molRegion.r_arr();
            const auto rOuter = molOuter.r_arr();
            const auto d2 = std::pow(rRegion[0] - rOuter[0], 2) + std::pow(rRegion[1] - rOuter[1], 2) + std::pow(rRegion[2] - rOuter[2], 2);
            if (d2 > cutoff2) continue;

            int bin = std::clamp((int) (std::sqrt(d2) / binSize), 0, numBins-1);
            threadCounts[mardyn_get_thread_num()][bin]++;
        }
        countedParticles++;
    }

    for (auto& buffer : threadCounts) {
        for (int idx = 0; idx < numBins; idx++) counts[idx] += buffer[idx];
    }
}

std::vector<double> TheoreticalRDF::getCurrentRDF() {
    if (countedParticles == 0) return std::vector<double>(numBins, 0.0);

    std::vector<double> result (numBins, 0.0);
    const auto factor = countedParticles * rhoGlobal;

    #pragma omp parallel for
    for (int idx = 0; idx < numBins; idx++) {
        const double rMin = (idx + 0) * binSize;
        const double rMax = (idx + 1) * binSize;
        const double binVol = (4.0/3.0) * M_PI * (std::pow(rMax, 3) - std::pow(rMin, 3));

        result[idx] = counts[idx] / (binVol * factor);
    }

    return result;
}

void TheoreticalRDF::writeFile(const std::string &name, unsigned long simstep) {
    if (simstep%1000 != 0) return;

    std::string fileName = name + std::to_string(simstep) + ".txt";
    std::ofstream outputFile(fileName);

    const auto rdf = getCurrentRDF();
    for (int idx = 0; idx < numBins; idx++) {
        outputFile << idx * binSize + binSize/2 << "\t" << rdf[idx] << std::endl;
    }

    outputFile.flush();
    outputFile.close();
}
