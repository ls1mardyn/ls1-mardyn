#ifndef DECOMPWRITER_H_
#define DECOMPWRITER_H_

#include "io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <fstream>
#include <sstream>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
//class Molecule;
//class DecompWriter; 

//! @todo comment
class DecompWriter : public OutputBase {
public:
	DecompWriter(unsigned long writeFrequency, std::string mode, std::string filename, unsigned long numberOfTimesteps, bool incremental);
	~DecompWriter();
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
	unsigned long _numberOfTimesteps;
	unsigned long _writeFrequency;
	std::string _mode;
	bool _filenameisdate;
	bool _incremental;
	std::string _filename;
};

#endif /*DECOMPWRITER_H_*/
