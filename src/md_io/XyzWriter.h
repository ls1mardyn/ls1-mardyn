#ifndef XYZWRITER_H_
#define XYZWRITER_H_

#include "md_io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>


class ParticleContainer;
class DomainDecompBase; 
class Domain;


using namespace std;

//! @todo comment
class XyzWriter : public OutputBase{
 public:
  XyzWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental);
  ~XyzWriter();
  //! @todo comment
  void initOutput(ParticleContainer* particleContainer,
		  DomainDecompBase* domainDecomp, Domain* domain);
  //! @todo comment
  void doOutput(
     ParticleContainer* particleContainer,
     DomainDecompBase* domainDecomp, Domain* domain,
     unsigned long simstep, list<ChemicalPotential>* lmu
  );
  //! @todo comment  
  void finishOutput(ParticleContainer* particleContainer,
		    DomainDecompBase* domainDecomp, Domain* domain);
 private:
  unsigned long _numberOfTimesteps; 
  unsigned long _writeFrequency;
  bool _filenameisdate;
  bool _incremental;
  string _filename;
};

#endif /*XYZWRITER_H_*/
