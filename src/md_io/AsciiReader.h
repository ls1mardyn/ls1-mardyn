#ifndef ASCIIREADER_H_
#define ASCIIREADER_H_

#include "md_io/InputBase.h"
#include <string>
#include <sstream>
#include <fstream>

class ChemicalPotential;

//class ParticleContainer;
//class DomainDecompBase;
//class Domain;

using namespace std;

class AsciiReader : public InputBase{
 public:
  AsciiReader();
  ~AsciiReader();
  void setPhaseSpaceFile(string filename);
  void setPhaseSpaceHeaderFile(string filename);
  void readPhaseSpaceHeader(Domain* domain, double timestep);
  unsigned long readPhaseSpace(ParticleContainer* particleContainer, list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);
 private:
#ifdef PARALLEL
  istringstream _phaseSpaceFileStream;
#else
  fstream _phaseSpaceFileStream;
#endif
  fstream _phaseSpaceHeaderFileStream;
  string _phaseSpaceFileName;
  string _phaseSpaceHeaderFileName;
};

#endif /*ASCIIREADER_H_*/
