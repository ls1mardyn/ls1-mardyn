#ifndef XYZWRITER_H_
#define XYZWRITER_H_

#include "io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

//! @todo comment
class XyzWriter : public OutputBase {
public:
	XyzWriter(unsigned long writeFrequency, std::string filename, unsigned long numberOfTimesteps, bool incremental);
	~XyzWriter();
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
	std::string _filename;
	unsigned long _numberOfTimesteps;
	unsigned long _writeFrequency;
	bool _filenameisdate;
	bool _incremental;
};

#endif /*XYZWRITER_H_*/
