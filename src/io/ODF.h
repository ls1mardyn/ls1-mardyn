/**
 * @file ODF.h
 * @author Joshua Marx
 * @date 09/10/19
 */

#pragma once

#include <particleContainer/adapter/ODFCellProcessor.h>
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"

/**
 * Calculates the orientation distribution function (ODF) for pure fluids or binary mixtures
 * Orientation distribution function is calculated in regards to the centers of the molecules,
 * site-wise ODFs are not implemented
 * This plugin is designed for use with Stockmayer fluids (LJ center + dipole), mixtures of Stockmayer
 * fluids with LJ fluids and Stockmayer fluid mixtures. Application to other rotating fluids is theoretically
 * possible, but requires some minor adjustments and more testing.
 *
 * For now, the plugin only works for molecules with exactly one dipole
 * @note IMPORTANT: the standard dipole orientation in the input must be set to [0 0 1]!
 * @note IMPORTANT: a molecule consisting of a LJ center and a dipole site in the same location will NOT rotate. Add the site 'Stockmayer' to the molecule to enable rotation with Ixx = Iyy = Izz = 1
 */
class ODF : public PluginBase {
public:
	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
	void readXML(XMLfileUnits& xmlconfig) override;
	void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
					 unsigned long simstep) override;
	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override;
	void reset();
	void collect(DomainDecompBase* domainDecomp);
	void calculateOrientation(const std::array<double, 3> &simBoxSize,
                              const Molecule &mol1,
                              const Molecule &mol2,
                              const std::array<double, 3> &orientationVector1);
	void output(Domain* domain, long unsigned timestep);
	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{};
	std::string getPluginName() override { return std::string("ODF"); }
	static PluginBase* createInstance() { return new ODF(); }

private:
	bool _readStartingStep;
	unsigned long _writeFrequency;      // Write frequency for all profiles -> Length of recording frame before output
	unsigned long _initStatistics;      // Timesteps to skip at start of the simulation
	unsigned long _recordingTimesteps;  // Record every Nth timestep during recording frame
	unsigned _phi1Increments;
	unsigned _phi2Increments;
	unsigned _gammaIncrements;
	unsigned _numPairs;
	unsigned _numComponents;
	unsigned long _numElements;
	double _shellCutOff;

	std::string _outputPrefix;

	std::vector<unsigned long> _ODF11;
	std::vector<unsigned long> _ODF12;
	std::vector<unsigned long> _ODF21;
	std::vector<unsigned long> _ODF22;
	// thread buffers to avoid race conditions
	std::vector<std::vector<unsigned long>> _threadLocalODF11;
	std::vector<std::vector<unsigned long>> _threadLocalODF12;
	std::vector<std::vector<unsigned long>> _threadLocalODF21;
	std::vector<std::vector<unsigned long>> _threadLocalODF22;

	std::unique_ptr<ODFCellProcessor> _cellProcessor;
};
