// Created by Joshua Marx 08/2019
// Only works for pure fluids and binary mixtures for now

#ifndef SRC_PLUGINS_ODF_H_
#define SRC_PLUGINS_ODF_H_

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"

class ODF : public PluginBase {
public:
	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
	void readXML(XMLfileUnits& xmlconfig) override;
	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override;
	void reset();
	void record(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp,
				unsigned long simstep);
	void calculateOrientation(ParticleContainer* particleContainer, Domain* domain, const Molecule& mol1);
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
	unsigned long _numComp1;
	unsigned long _numComp2;
	unsigned long _numElements;
	unsigned long _numMolecules;
	double _shellCutOff[2];
	bool _mixingRule;
	std::string _outputPrefix;
//	std::map<unsigned long, unsigned> _cid;
//	std::map<unsigned long, unsigned> _molID;
//	std::map<unsigned long, unsigned long> _doODF;
//	std::map<unsigned long, std::map<unsigned long, double>> _Q1;
//	std::map<unsigned long, std::map<unsigned long, double>> _upVec1;
//	std::map<unsigned long, std::map<unsigned long, double>> _r1;
//	std::map<unsigned long, unsigned long> _ODF11;
//	std::map<unsigned long, unsigned long> _ODF12;
//	std::map<unsigned long, unsigned long> _ODF21;
//	std::map<unsigned long, unsigned long> _ODF22;
	std::map<unsigned long, unsigned long> _localODF11;
	std::map<unsigned long, unsigned long> _localODF12;
	std::map<unsigned long, unsigned long> _localODF21;
	std::map<unsigned long, unsigned long> _localODF22;
};

#endif /* SRC_PLUGINS_ODF_H_ */
