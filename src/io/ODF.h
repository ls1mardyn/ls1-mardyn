// Created by Joshua Marx 08/2019
// Only works for pure fluids and binary mixtures for now

#pragma once

#include <particleContainer/adapter/ODFCellProcessor.h>
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"

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
	void calculateOrientation(const array<double, 3> &simBoxSize,
                              const Molecule &mol1,
                              const Molecule &mol2,
                              const array<double, 3> &upVec1);
	void output(Domain* domain, long unsigned timestep);
	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{};
	std::string getPluginName() override { return std::string("ODF"); }
	static PluginBase* createInstance() { return new ODF(); }

private:
	unsigned long _writeFrequency;      // Write frequency for all profiles -> Length of recording frame before output
	unsigned long _initStatistics;      // Timesteps to skip at start of the simulation
	unsigned long _recordingTimesteps;  // Record every Nth timestep during recording frame
	unsigned _phi1Increments;
	unsigned _phi2Increments;
	unsigned _gammaIncrements;
	unsigned _numPairs;
	unsigned _numComponents;
	unsigned long _numElements;
	double _shellCutOff[2];
	bool _mixingRule;
	std::string _outputPrefix;

	std::vector<unsigned long> _ODF11;
	std::vector<unsigned long> _ODF12;
	std::vector<unsigned long> _ODF21;
	std::vector<unsigned long> _ODF22;
	std::vector<unsigned long> _localODF11;
	std::vector<unsigned long> _localODF12;
	std::vector<unsigned long> _localODF21;
	std::vector<unsigned long> _localODF22;

	std::unique_ptr<ODFCellProcessor> _cellProcessor;
};
