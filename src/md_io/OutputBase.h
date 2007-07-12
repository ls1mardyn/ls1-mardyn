#ifndef OUTPUTBASE_H_
#define OUTPUTBASE_H_

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
  class OutputBase;
}
using namespace std;

//! @brief interface for any kind of output class
//!
//! @todo more comment
class md_io::OutputBase{
 public:
  OutputBase(){}
  
  virtual ~OutputBase(){}
  
  //! @brief to be called at the beginning of the simulation
  //!
  //! @todo comment 
  virtual void initOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain) = 0;
  
  //! @brief to be called in each time step
  //!
  //! @todo comment
  virtual void doOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) = 0;
  
  //! @brief to be called at the end of the simulation
  //!
  //! @todo comment
  virtual void finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain) = 0;
  
};

#endif /*OUTPUTBASE_H_*/
