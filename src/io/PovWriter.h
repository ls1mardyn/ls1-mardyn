#ifndef POVWRITER_H_
#define POVWRITER_H_

#include "io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

class PovWriter : public OutputBase {
public:
	//! @brief writes a POVray file of the current state of the simluation
	//!
	//! @param filename Name of the POV file (including path)
	PovWriter(unsigned long writeFrequency, std::string filename, unsigned long numberOfTimesteps, bool incremental);
	~PovWriter();
	void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	void doOutput(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain,
			unsigned long simstep, std::list<ChemicalPotential>* lmu
	);
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
private:
	std::string _filename;
	unsigned long _writeFrequency;
	unsigned long _numberOfTimesteps;
	bool  _incremental;
	bool  _filenameisdate;
};

#endif /*POVWRITER_H_*/
