/*
 * FlopRateWriter.h
 *
 *  Created on: 8 May 2017
 *      Author: tchipevn
 */

#ifndef SRC_IO_FLOPRATEWRITER_H_
#define SRC_IO_FLOPRATEWRITER_H_

#include "OutputBase.h"
class FlopCounter;

class FlopRateWriter: public OutputBase {
public:
	FlopRateWriter(double cutoff, double LJcutoff) :
			OutputBase(), _cutoffRadius(cutoff), _LJCutoffRadius(LJcutoff),
			_writeToStdout(false), _writeToFile(false) {
	}
	~FlopRateWriter() {}

	//! @brief will be called at the beginning of the simulation
	//!
	//! Some OutputPlugins will need some initial things to be done before
	//! the output can start, e.g. opening some files. This method will
	//! be called once at the beginning of the simulation (see Simulation.cpp)
	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

	void readXML(XMLfileUnits& /*xmlconfig*/);

	//! @brief will be called in each time step
	//!
	//! Most of the times, the output should either be done every time step or at least
	//! every n-th time step. Therefore, this method has an additional parameter simstep,
	//! allowing to do a output depending on the current simulation time step. This method
	//! will be called once every time step during the simulation (see Simulation.cpp)
	void doOutput(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep,
			std::list<ChemicalPotential>* lmu,
			std::map<unsigned, CavityEnsemble>* mcav
	);

	//! @brief will be called at the end of the simulation
	//!
	//! Some OutputPlugins will need to do some things at the end of the simulation,
	//! e.g. closing some files. This method will
	//! be called once at the end of the simulation (see Simulation.cpp)
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

	std::string getPluginName() {
		return std::string("FlopRateOutputPlugin.");
	}

private:
	double _cutoffRadius, _LJCutoffRadius;

	bool _writeToStdout, _writeToFile;
	std::ofstream _fileStream;
	unsigned long _writeFrequency;
	std::string _outputPrefix;
};

#endif /* SRC_IO_FLOPRATEWRITER_H_ */
