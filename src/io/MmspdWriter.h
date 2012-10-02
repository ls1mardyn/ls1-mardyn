#ifndef MMSPDWRITER_H_
#define MMSPDWRITER_H_

// by Stefan Becker <stefan.becker@mv.uni-kl.de>
//output writer wirting the output in the *.mmspd file format required for MegaMol
// a detailed documentation of this file format can be obtained from the MegaMol documentation (see website)

#include "io/OutputBase.h"
#include "ensemble/GrandCanonical.h"
#include <string>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

class MmspdWriter : public OutputBase{
  //!@brief: writes a mmspd file used by MegaMol
 
  public:
    MmspdWriter(unsigned long writeFrequency, std::string filename, unsigned long numberOfTimesteps, bool incremental);
	~MmspdWriter();
    
    void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);
	//! @todo comment
    void doOutput( ParticleContainer* particleContainer,
		   DomainDecompBase* domainDecomp, Domain* domain,
		   unsigned long simstep, std::list<ChemicalPotential>* lmu
	);
	//! @todo comment
    void finishOutput( ParticleContainer* particleContainer,
		       DomainDecompBase* domainDecomp, Domain* domain);
  private:
      std::string _filename;
      unsigned long _numberOfTimesteps;
      unsigned long _writeFrequency;
      bool _filenameisdate;
      bool _incremental;
};

#endif /*MMSPDWRITER_H_*/
