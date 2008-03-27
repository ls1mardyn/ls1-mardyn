#ifndef INPUTBASE_H_
#define INPUTBASE_H_

#include "Domain.h"
#include <string>

namespace datastructures{
  template<class ParticleType>
  class ParticleContainer;
}

namespace parallel{
  class DomainDecompBase; 
}

class Molecule;

namespace md_io{
  class InputBase;
}
using namespace std;

//! @brief interface for any kind of input class
//!
//! @todo more comment
class md_io::InputBase{
 public:
  InputBase(){}
  
  virtual ~InputBase(){}
  
  //! @brief set the phase space file name
  virtual void setPhaseSpaceFile(string filename) = 0;
  
  //! @brief set the phase space header file name (can be identical to the
  //         phase space file
  virtual void setPhaseSpaceHeaderFile(string filename) = 0;
  
  //! @brief read the phase space components and header information
  virtual void readPhaseSpaceHeader(Domain* domain) = 0;
  
  //! @brief read the actual phase space information
  virtual void readPhaseSpace(datastructures::ParticleContainer<Molecule>* particleContainer, Domain* domain) = 0;
  
};

#endif /*INPUTBASE_H_*/
