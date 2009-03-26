#ifndef POVWRITER_H_
#define POVWRITER_H_

#include "md_io/OutputBase.h"
//#include "molecules/Component.h"
#include <string>
//#include <fstream>
//#include <sstream>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

using namespace std;

class PovWriter : public OutputBase{
 public:
   //! @brief writes a POVray file of the current state of the simluation
   //!
   //! @param filename Name of the POV file (including path)
  PovWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental);
  ~PovWriter();
  void initOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain);
  void doOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);
  void finishOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain);
 private:
  string _filename;
  unsigned long _writeFrequency;
  unsigned long _numberOfTimesteps;
  bool  _incremental;
  bool  _filenameisdate;
};

#endif /*POVWRITER_H_*/
