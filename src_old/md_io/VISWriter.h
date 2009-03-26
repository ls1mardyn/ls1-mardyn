#ifndef VISWRITER_H_
#define VISWRITER_H_

#include "md_io/OutputBase.h"
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
  class VISWriter; 
}
using namespace std;

class md_io::VISWriter : public md_io::OutputBase{
 public:
   //! @brief writes a VIS file
   //!
  VISWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental);
  ~VISWriter();
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
  bool	_incremental;
  bool	_filenameisdate;
};

#endif /*VISWRITER_H_*/
