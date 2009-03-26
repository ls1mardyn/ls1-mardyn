#ifndef POVWRITER_H_
#define POVWRITER_H_

#include "md_io/OutputBase.h"
#include "molecules/Component.h"
#include <string>
#include <fstream>

namespace datastructures{
  template<class ParticleType>
  class ParticleContainer;
}

namespace parallel{
  class DomainDecompBase; 
}

class Domain;
class Molecule;

namespace md_io{
  class PovWriter; 
}
using namespace std;

class md_io::PovWriter : public md_io::OutputBase{
 public:
   //! @brief writes a POVray file of the current state of the simluation
   //!
   //! @param filename Name of the POV file (including path)
  PovWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental);
  ~PovWriter();
  void initOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain);
  void doOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);
  void finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain);
 private:
  string _filename;
  unsigned long _writeFrequency;
  unsigned long _numberOfTimesteps;
  bool  _incremental;
  bool  _filenameisdate;
};

#endif /*POVWRITER_H_*/
