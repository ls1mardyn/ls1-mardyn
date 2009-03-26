#ifndef XMLREADER_H_
#define XMLREADER_H_

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
  class XMLReader; 
}
using namespace std;

class md_io::XMLReader : public md_io::InputBase{
 public:
  XMLReader();
  ~XMLReader();
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

#endif /*XMLREADER_H_*/
