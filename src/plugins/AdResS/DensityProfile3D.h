//
// Created by alex on 29.12.23.
//

#ifndef MARDYN_DENSITYPROFILE3D_H
#define MARDYN_DENSITYPROFILE3D_H

#include <array>
#include <vector>
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

class DensityProfile3D {
public:
    void init(double binWidth, Domain* domain);
    void sampleDensities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
    [[nodiscard]] const std::vector<double>& getDensity(int dim) const;
private:
    double _binWidth;
    std::array<unsigned long, 3> _binDims;
    std::array<double, 3> _binVolumes;
    std::array<std::vector<double>, 3> _localDensities;
    std::array<std::vector<double>, 3> _globalDensities;
    void resetBuffers();

    void sampleForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
};


#endif //MARDYN_DENSITYPROFILE3D_H
