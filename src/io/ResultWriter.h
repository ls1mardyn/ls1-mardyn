#ifndef RESULTWRITER_H_
#define RESULTWRITER_H_

#include <fstream>
#include <string>

#include "ensemble/GrandCanonical.h"
#include "io/OutputBase.h"
#include "utils/Accumulator.h"


//! @brief writes thermodynamic properties to a file
//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
//!
//! The following global values will be written to a file:
//! - Simulation time step
//! - time since the simulation started (dimensionless)
//! - Average potential Energy
//! - Pressure
//! - BetaTrans
//! - BetaRot
//!
//! @param writeFrequency Controls the frequency of writing out the data
//! (every timestep, every 10th, 100th, ... timestep)
//!
class ResultWriter : public OutputBase {
public:
	ResultWriter(){}
	ResultWriter(unsigned long writeFrequency, std::string outputPrefix);
	~ResultWriter();

	virtual void readXML(XMLfileUnits& xmlconfig);

	//! @todo comment
	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	//! @todo comment
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu
	);
	//! @todo comment
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("ResultWriter");
	}

private:
	//! prefix for the names of all output files
	std::ofstream _resultStream;
	unsigned long _writeFrequency;
	std::string _outputPrefix;
	Accumulator<double> *_U_pot_acc;
	Accumulator<double> *_p_acc;
};

#endif /*RESULTWRITER_H_*/
