#ifndef CHECKPOINTWRITER_H_
#define CHECKPOINTWRITER_H_

#include "md_io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>

class ParticleContainer;
class DomainDecompBase; 
class Domain;




using namespace std;

class CheckpointWriter : public OutputBase{
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
  void initOutput(ParticleContainer* particleContainer,
		  DomainDecompBase* domainDecomp, Domain* domain);
  void doOutput(
     ParticleContainer* particleContainer,
     DomainDecompBase* domainDecomp, Domain* domain,
     unsigned long simstep, list<ChemicalPotential>* lmu
  );
  void finishOutput(ParticleContainer* particleContainer,
		    DomainDecompBase* domainDecomp, Domain* domain);
 private:
  string _filename;
  unsigned long _writeFrequency;
  unsigned long _numberOfTimesteps;
  bool	_incremental;
  bool	_filenameisdate;
};

#endif /*CHECKPOINTWRITER_H_*/
