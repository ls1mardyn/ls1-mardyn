#ifndef CHECKPOINTWRITER_H_
#define CHECKPOINTWRITER_H_

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
  class CheckpointWriter; 
}
using namespace std;

class md_io::CheckpointWriter : public md_io::OutputBase{
 public:
   //! @brief writes a checkpoint file that can be used to continue the simulation
   //!
   //! The format of the checkpointfile written by this method is the same as the format
   //! of the input file.
   //! @param filename Name of the checkpointfile (including path)
   //! @param particleContainer The molecules that have to be written to the file are stored here
   //! @param domainDecomp In the parallel version, the file has to be written by more than one process.
   //!                     Methods to achieve this are available in domainDecomp
  CheckpointWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental);
  ~CheckpointWriter();
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

#endif /*CHECKPOINTWRITER_H_*/
