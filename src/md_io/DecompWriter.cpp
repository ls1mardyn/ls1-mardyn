// DecompWriter.cpp

#include "md_io/DecompWriter.h"
#include "Common.h"
#include "Domain.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"

DecompWriter::DecompWriter(unsigned long writeFrequency, string mode, string filename, unsigned long numberOfTimesteps, bool incremental) {
  _filename = filename;
  _mode = mode;
  _writeFrequency = writeFrequency;
  _incremental = incremental;
  _numberOfTimesteps = numberOfTimesteps;

  if (filename == "default")
    _filenameisdate = true;
  else
    _filenameisdate = false;
}

DecompWriter::~DecompWriter(){}

void DecompWriter::initOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
}

void DecompWriter::doOutput( ParticleContainer* particleContainer,
                             DomainDecompBase* domainDecomp, Domain* domain,
			     unsigned long simstep, list<ChemicalPotential>* lmu ) 
{
  if(simstep%_writeFrequency == 0) {
    stringstream filenamestream;
    if(_filenameisdate) {
      filenamestream << gettimestring() << ".out";
    } else {
      filenamestream << _filename;
    }

    if(_incremental) {
      unsigned long temp = simstep/_writeFrequency;
      filenamestream << "-";
      while(temp < floor((double) (_numberOfTimesteps/_writeFrequency))){
        filenamestream << "0";
        temp = temp*10;
      }
      filenamestream << simstep/_writeFrequency << ".decomp";
    } else {
      filenamestream << ".decomp";
    }

    domainDecomp->printDecomp(filenamestream.str(), domain);
    
    if(_mode=="withParticles"){
      int ownRank = domainDecomp->getRank();
      for(int process = 0; process < domainDecomp->getNumProcs(); process++){
        if(ownRank==process){
          ofstream decompstrm(filenamestream.str().c_str(), ios::app);
          if(ownRank==0) decompstrm << "particleData xyz" << endl;
          Molecule* tempMol;
          for(tempMol = particleContainer->begin(); tempMol != particleContainer->end(); tempMol = particleContainer->next()){
            decompstrm << tempMol->r(0) << "\t" << tempMol->r(1) << "\t" << tempMol->r(2) << endl;
          }
          decompstrm.close();
        }
        domainDecomp->barrier();
      }
    }
    else if(domainDecomp->getRank()==0){
      ofstream decompstrm(filenamestream.str().c_str(), ios::app);
      decompstrm << "particleData none" << endl;
      decompstrm.close();
    }
  }  
}

void DecompWriter::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
}
