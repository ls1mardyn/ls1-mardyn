//
// Created by alex on 29.12.23.
//

#ifndef MARDYN_DENSITYPROFILE3D_H
#define MARDYN_DENSITYPROFILE3D_H

#include <array>
#include <vector>
#include <string>
#include "Domain.h"
#include "Interpolation.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

class DensityProfile3D {
public:
	enum Type {
		SAMPLE, SMOOTH, GMM
	};
    void init(double binWidth, Domain* domain, double smoothingFactor);
    void sampleDensities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
    [[nodiscard]] const std::vector<double>& getDensity(int dim) const;
    [[nodiscard]] std::vector<double> getDensitySmoothed(int dim) const;
    void computeGMMDensities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
    [[nodiscard]] const Interpolation::Function& getGMMDensity(int dim) const;
	//! @brief writes densities currently in global buffer
	void writeDensity(const std::string &filename, const std::string &separator, int dim, Type type);
private:
	double _smoothingFactor;
    double _binWidth;
    std::array<unsigned long, 3> _binDims;
    std::array<double, 3> _binVolumes;
    std::array<std::vector<double>, 3> _localDensities;
    std::array<std::vector<double>, 3> _globalDensities;
    std::array<Interpolation::Function, 3> _gmmDensities;
    Interpolation::Matrix _smoothingFilter;
    void resetBuffers();
};


#endif //MARDYN_DENSITYPROFILE3D_H
