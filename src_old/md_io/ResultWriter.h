#ifndef RESULTWRITER_H_
#define RESULTWRITER_H_

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
  class ResultWriter; 
}
using namespace std;

//! @brief writes thermodynamic properties to a file
//! @author Martin Buchholz
//!
//! The following global values will be written to a file:
//! - Simulation time step
//! - time since the simulation started (dimensionless)
//! - Average potential Energy
//! - Pressure
//! - BetaTrans
//! - BetaRot
class md_io::ResultWriter : public md_io::OutputBase{
 public:
  ResultWriter(string outputPrefix);
  ~ResultWriter();
  //! @todo comment
  void initOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain);
  //! @todo comment
  void doOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);
  //! @todo comment
  void finishOutput(datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp, Domain* domain);
  
 private:
  //! prefix for the names of all output files
  string _outputPrefix;
  ofstream _resultStream;
};

#endif /*RESULTWRITER_H_*/
