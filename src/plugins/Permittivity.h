/*
 * Permittivity.h
 *
 *  Created on: November 2019
 *      Author: Joshua Marx
 */
 
 // DESCRIPTION: Samples the relative permittivity of Stockmayer fluids in the NVT ensemble
 // Important: Permittivity may take a long time to converge, i.e. a few million steps with ~1000 particles. Reducing number of slabs for the thermostat can drastically improve results!
 
#ifndef SRC_PLUGINS_PERMITTIVITY_H_
#define SRC_PLUGINS_PERMITTIVITY_H_

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"

class Permittivity: public PluginBase {
public:

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
	void readXML(XMLfileUnits& xmlconfig) override;
	void record(ParticleContainer* particleContainer);
	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override;
	void reset();
	void collect(DomainDecompBase* domainDecomp);
	void writeRunningAverage(unsigned long indexM, double tempMX, double tempMY, double tempMZ, double tempMSquared);
	void output(Domain* domain, unsigned long timestep);
	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{};
	std::string getPluginName() override { return std::string("Permittivity"); }
	static PluginBase* createInstance() { return new Permittivity();}
	
	
	
private:
	bool _readStartingStep;             // Auxiliary bool variable to read the current time step during the first iteration of endStep
	unsigned long _writeFrequency;      // Write frequency for all profiles -> Length of recording frame before output
	unsigned long _initStatistics;      // Timesteps to skip at start of the simulation
	unsigned long _recordingTimesteps;  // Record every Nth timestep during recording frame
	unsigned long _accumulatedSteps;    // Number of steps in a block
	unsigned long _runningAverageSteps; // Timesteps for running average output
	unsigned long _totalNumTimeSteps;	// Total number of time steps
	unsigned long _startingStep;        // current time step at the start of simulation
	unsigned long _numParticlesLocal;   // Counts number of particles considered for each block locally
	unsigned long _RAVCounter;
	unsigned _numOutputs;               // Number of output values to be written
	unsigned _numComponents;            // Number of components
	unsigned _currentOutputNum;         // Number of current block
	std::map<unsigned,double> _myAbs;   // Dipole moment of component
	std::map<unsigned,unsigned long> _numParticles; // Counts number of particles considered for each block globally
	std::map<unsigned,double> _outputSquaredM; // Total average squared dipole moment for each block
	std::map<unsigned, std::map<unsigned, double>> _outputM; // Total average dipole moment for each block
	std::map<unsigned long, std::map<unsigned long, double >> _localM; // Total dipole moment local
	std::map<unsigned long, std::map<unsigned long, double >> _globalM; // Total dipole moment global
	std::array<double, 3> simBoxSize;
	double _totalAverageM[3]; // Total dipole moment averaged over whole production run
	double _totalAverageSquaredM; // Total squared dipole moment averaged over whole production run
	double _MRAV[3]; // Running average for M
	double _MSquaredRAV; // Running average for <M2>
	double _permRAV; // Running average for permittivity (without consideration of <M>2
	double _V; // Box volume
	double _T; // Target temperature
	std::ofstream _ravStream;
	std::string _outputPrefix;
};

#endif /*SRC_PLUGINS_PERMITTIVITY_H_*/
