#ifndef INPUTBASE_H_
#define INPUTBASE_H_

#include <string>

class ParticleContainer;
class DomainDecompBase;
class Domain;

using namespace std;

//! @brief interface for any kind of input class
//!
//! @todo more comment
class InputBase{
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
  virtual void readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) = 0;
  
};

#endif /*INPUTBASE_H_*/
