#ifndef DECOMPWRITER_H_
#define DECOMPWRITER_H_

#include "md_io/OutputBase.h"
#include <string>
#include <fstream>
#include <sstream>

class ParticleContainer;
class DomainDecompBase; 
class Domain;
//class Molecule;
//class DecompWriter; 


using namespace std;

//! @todo comment
class DecompWriter : public OutputBase{
 public:
  DecompWriter(unsigned long writeFrequency, string mode, string filename, unsigned long numberOfTimesteps, bool incremental);
  ~DecompWriter();
  //! @todo comment
  void initOutput(ParticleContainer* particleContainer,
		  DomainDecompBase* domainDecomp, Domain* domain);
  //! @todo comment
  void doOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);
  //! @todo comment  
  void finishOutput(ParticleContainer* particleContainer,
		    DomainDecompBase* domainDecomp, Domain* domain);
 private:
  unsigned long _numberOfTimesteps; 
  unsigned long _writeFrequency;
  string _mode;
  bool _filenameisdate;
  bool _incremental;
  string _filename;
};

#endif /*DECOMPWRITER_H_*/
