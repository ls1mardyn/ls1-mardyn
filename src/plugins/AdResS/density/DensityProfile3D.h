//
// Created by alex on 29.12.23.
//

#ifndef MARDYN_DENSITYPROFILE3D_H
#define MARDYN_DENSITYPROFILE3D_H

#include <array>
#include <vector>
#include <string>
#include "Domain.h"
#include "plugins/AdResS/Interpolation.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "GridSampler.h"
#include "plugins/AdResS/util/Averager.h"

class DensityProfile3D {
public:
	enum Type {
		SAMPLE, SMOOTH, GMM, FT, GRID
	};
    void init(double binWidth, Domain* domain, double smoothingFactor, double gridRadius);
    void sampleDensities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
    [[nodiscard]] const std::vector<double>& getDensity(int dim) const;
    [[nodiscard]] std::vector<double> getDensitySmoothed(int dim) const;
    void computeGMMDensities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
	void computeFTDensities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
    [[nodiscard]] const Interpolation::Function& getGMMDensity(int dim) const;
    [[nodiscard]] const Interpolation::Function& getFTDensity(int dim) const;
	//! @brief writes densities currently in global buffer
	void writeDensity(const std::string &filename, const std::string &separator, int dim, Type type);
	void step(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
private:
	double _smoothingFactor;
    double _binWidth;
    std::array<unsigned long, 3> _binDims;
    std::array<double, 3> _binVolumes;
    std::array<std::vector<double>, 3> _localDensities;
    std::array<std::vector<double>, 3> _globalDensities;
	std::array<Interpolation::Function, 3> _ftDensities;
    std::array<Interpolation::Function, 3> _gmmDensities;
    Interpolation::Matrix _smoothingFilter;
    void resetBuffers();
	std::array<std::vector<double>, 3> getGlobalMolPos(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
	Grid3D::GridSampler _gridSampler;
	Grid3D::Grid _grid;
	Averager<std::vector<double>> _averager;
};


#endif //MARDYN_DENSITYPROFILE3D_H
