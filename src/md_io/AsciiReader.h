#ifndef ASCIIREADER_H_
#define ASCIIREADER_H_

#include "md_io/InputBase.h"
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
  class AsciiReader; 
}
using namespace std;

class md_io::AsciiReader : public md_io::InputBase{
 public:
  AsciiReader();
  ~AsciiReader();
  void setPhaseSpaceFile(string filename);
  void setPhaseSpaceHeaderFile(string filename);
  void readPhaseSpaceHeader(Domain* domain);
  void readPhaseSpace(datastructures::ParticleContainer<Molecule>* particleContainer, Domain* domain);
 private:
  fstream _phaseSpaceFileStream;
  fstream _phaseSpaceHeaderFileStream;
  string _phaseSpaceFileName;
  string _phaseSpaceHeaderFileName;
};

#endif /*ASCIIREADER_H_*/
