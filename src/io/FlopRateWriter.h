/*
 * FlopRateWriter.h
 *
 *  Created on: 8 May 2017
 *      Author: tchipevn
 */

#ifndef SRC_IO_FLOPRATEWRITER_H_
#define SRC_IO_FLOPRATEWRITER_H_

#include "plugins/PluginBase.h"

#include <fstream>
#include <sstream>

class FlopCounter;

class FlopRateWriter: public PluginBase {
public:
	FlopRateWriter() : PluginBase(), _flopCounter(nullptr) {}
	~FlopRateWriter() {}

	//! @brief will be called at the beginning of the simulation
	//!
	//! Some OutputPlugins will need some initial things to be done before
	//! the output can start, e.g. opening some files. This method will
	//! be called once at the beginning of the simulation (see Simulation.cpp)
	void init(ParticleContainer *particleContainer,
              DomainDecompBase *domainDecomp, Domain *domain);

	void readXML(XMLfileUnits& /*xmlconfig*/);

	//! @brief will be called after forces have been applied and exchanged

	void afterForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep);

	//! @brief will be called in each time step
	//!
	//! Most of the times, the output should either be done every time step or at least
	//! every n-th time step. Therefore, this method has an additional parameter simstep,
	//! allowing to do a output depending on the current simulation time step. This method
	//! will be called once every time step during the simulation (see Simulation.cpp)
	void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep
    );

	//! @brief will be called at the end of the simulation
	//!
	//! Some OutputPlugins will need to do some things at the end of the simulation,
	//! e.g. closing some files. This method will
	//! be called once at the end of the simulation (see Simulation.cpp)
	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() {
		return std::string("FlopRateWriter");
	}
 	static PluginBase* createInstance() { return new FlopRateWriter(); }

	void measureFLOPS(ParticleContainer* particleContainer, unsigned long simstep);

private:
	void setPrefix(double f_in, double& f_out, char& prefix) const;

	FlopCounter * _flopCounter;
	bool _writeToStdout, _writeToFile;
	std::ofstream _fileStream;
	unsigned long _writeFrequency;
	std::string _outputPrefix;
};

#endif /* SRC_IO_FLOPRATEWRITER_H_ */
