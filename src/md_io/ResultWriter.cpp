#include <ctime>
#include "md_io/ResultWriter.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"

ResultWriter::ResultWriter(string outputPrefix){
  _outputPrefix = outputPrefix;
}

ResultWriter::~ResultWriter(){}

void ResultWriter::initOutput(ParticleContainer* particleContainer,
			      DomainDecompBase* domainDecomp, Domain* domain){
   
  // initialize result file
  string resultfile(_outputPrefix+".res");
  time_t now;
  time(&now);
  if(domainDecomp->getRank()==0){
    _resultStream.open(resultfile.c_str());
    _resultStream << "# moldy MD simulation starting at " << ctime(&now) << endl;
    _resultStream << "#\tt\tU_pot\tPressure\tbeta_trans\tbeta_rot\t\tclock()/step" << endl;
    
  }
}

void ResultWriter::doOutput(ParticleContainer* particleContainer,
			    DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep){
  if(domainDecomp->getRank()==0){
    _resultStream << simstep << "\t" << domain->getCurrentTime()
                  << "\t" << domain->getAverageGlobalUpot() << "\t" << domain->getGlobalPressure()
                  << "\t" << domain->getGlobalBetaTrans() << "\t" << domain->getGlobalBetaRot() << endl;
  }
}

void ResultWriter::finishOutput(ParticleContainer* particleContainer,
				DomainDecompBase* domainDecomp, Domain* domain){
  _resultStream.close();
}
