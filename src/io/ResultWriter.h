#ifndef RESULTWRITER_H_
#define RESULTWRITER_H_

#include "io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <fstream>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

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
class ResultWriter : public OutputBase {
public:
	ResultWriter(unsigned long writeFrequency, std::string outputPrefix);
	~ResultWriter();
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

private:
	//! prefix for the names of all output files
	std::ofstream _resultStream;
	unsigned long _writeFrequency;
	std::string _outputPrefix;
};

#endif /*RESULTWRITER_H_*/
