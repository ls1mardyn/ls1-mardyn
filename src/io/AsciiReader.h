#ifndef ASCIIREADER_H_
#define ASCIIREADER_H_

#include <string>
#include <sstream>
#include <fstream>

#include "io/InputBase.h"


class ParticleContainer;
class DomainDecompBase;
class Domain;
class ChemicalPotential;

//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
class AsciiReader : public InputBase {
public:
	AsciiReader();
	~AsciiReader();
	void setPhaseSpaceFile(std::string filename);
	void setPhaseSpaceHeaderFile(std::string filename);
	void readPhaseSpaceHeader(Domain* domain, double timestep);
	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);
private:
#ifdef ENABLE_MPI
	std::istringstream _phaseSpaceFileStream;
#else
	std::fstream _phaseSpaceFileStream;
#endif
	std::fstream _phaseSpaceHeaderFileStream;
	std::string _phaseSpaceFileName;
	std::string _phaseSpaceHeaderFileName;
};

#endif /*ASCIIREADER_H_*/
