#ifndef GAMMAWRITER_H_
#define GAMMAWRITER_H_

#include "io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <fstream>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

//! @brief writes Surface tension to a file

class GammaWriter : public OutputBase {
public:
	GammaWriter(unsigned long writeFrequency, std::string outputPrefix);
	~GammaWriter();
	//! @todo comment
	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	//! @todo comment
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu,
			std::map<unsigned, CavityEnsemble>* mcav
	);
	//! @todo comment
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	
	std::string getPluginName() {
		return std::string("GammaWriter");
	}

private:
	//! prefix for the names of all output files
	std::ofstream _gammaStream;
	unsigned long _writeFrequency;
	std::string _outputPrefix;
};

#endif /*GAMMAWRITER_H_*/
