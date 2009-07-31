#ifndef VISWRITER_H_
#define VISWRITER_H_

#include "md_io/OutputBase.h"
#include "Domain.h"
#include "ensemble/GrandCanonical.h"
#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
class ChemicalPotential;

using namespace std;

class VISWriter : public OutputBase{
 public:
   //! @brief writes a .vis_ file
   //!
  VISWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental);
  ~VISWriter();
  void initOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain);
  void doOutput(
     ParticleContainer* particleContainer,
     DomainDecompBase* domainDecomp, Domain* domain,
     unsigned long simstep, list<ChemicalPotential>* lmu
  );
  void finishOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain);
 private:
  string _filename;
  unsigned long _writeFrequency;
  unsigned long _numberOfTimesteps;
  bool	_incremental;
  bool	_filenameisdate;
  bool  _wroteVIS;
};

#endif /*VISWRITER_H_*/
