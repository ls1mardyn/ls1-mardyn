#include "md_io/CheckpointWriter.h"
#include "md_io/Common.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include <sstream>


CheckpointWriter::CheckpointWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental) {
  _filename = filename;
  _writeFrequency = writeFrequency;
  _incremental = incremental;
  _numberOfTimesteps = numberOfTimesteps;

  if (filename == "default")
    _filenameisdate = true;
  else
    _filenameisdate = false;
}

CheckpointWriter::~CheckpointWriter(){}

void CheckpointWriter::initOutput(ParticleContainer* particleContainer,
          DomainDecompBase* domainDecomp, Domain* domain) {
}

void CheckpointWriter::doOutput(ParticleContainer* particleContainer,
        DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
  if(simstep%_writeFrequency == 0) {
    stringstream filenamestream;
    if(_filenameisdate) {
      filenamestream << gettimestring();
    } else {
      filenamestream << _filename;
    }
    if(_incremental) {
      unsigned long temp = simstep/_writeFrequency;
      filenamestream << "-";
      while(temp < floor(_numberOfTimesteps/_writeFrequency)){
  filenamestream << "0";
  temp = temp*10;
      }
      filenamestream << simstep/_writeFrequency << ".inp";
    } else {
      filenamestream << ".inp";
    }

    ofstream checkpointfilestream(filenamestream.str().c_str());
    checkpointfilestream << "MDProject\t20070111"<< endl;
    //! @todo use real cuurenttime
    checkpointfilestream << "currentTime\t0"  << endl;
    checkpointfilestream << "Temperature\t" << domain->getGlobalTemperature() << endl;
    checkpointfilestream << "Length\t" << domain->getGlobalLength(0) << "\t" 
                                       << domain->getGlobalLength(1) << "\t" 
                                       << domain->getGlobalLength(2) << endl;
    vector<Component> dcomponents = domain->getComponents();
    checkpointfilestream << "NumberOfComponents\t" << dcomponents.size() << endl;
    for(vector<Component>::const_iterator pos=dcomponents.begin(); pos!=dcomponents.end();++pos){
      pos->write(checkpointfilestream);
    }
    unsigned int numperline=dcomponents.size();
    unsigned int iout=0;
    vector<double> dmixcoeff = domain->getmixcoeff();
    for(vector<double>::const_iterator pos=dmixcoeff.begin(); pos!=dmixcoeff.end();++pos){
      checkpointfilestream << *pos;
      iout++;
      // 2 parameters (xi and eta)
      if(iout/2>=numperline) {
        checkpointfilestream << endl;
        iout=0;
        --numperline;
      }
      else if(!(iout%2)) {
        checkpointfilestream << "\t";
      }
      else {
        checkpointfilestream << " ";
      }
    }
    checkpointfilestream << domain->getepsilonRF() << endl;
    checkpointfilestream << "NumberOfMolecules\t" 
       << domain->getglobalNumMolecules() << endl;
      
    checkpointfilestream << "MoleculeFormat\t" << "\tICRVQD" << endl;
    checkpointfilestream.close();
      
    domainDecomp->writeMoleculesToFile(filenamestream.str().c_str(), particleContainer);

  }

}

void CheckpointWriter::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
}
