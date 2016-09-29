
#include "io/GammaWriter.h"

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"

#include <ctime>
using std::ios_base;

using namespace std;

GammaWriter::GammaWriter(unsigned long writeFrequency, string outputPrefix)
: _writeFrequency(writeFrequency),
	_outputPrefix(outputPrefix)
{ }

GammaWriter::~GammaWriter(){}

void GammaWriter::initOutput(ParticleContainer* /*particleContainer*/,
			      DomainDecompBase* domainDecomp, Domain* /*domain*/){
	 
	// initialize result file
	string resultfile(_outputPrefix+".gamma");
	_gammaStream.precision(6);
	time_t now;
	time(&now);
	if(domainDecomp->getRank()==0){
		_gammaStream.open(resultfile.c_str());
		_gammaStream << "# mardyn MD simulation starting at " << ctime(&now) << endl;
		_gammaStream << "#\tgamma" << endl;
	}
}

void GammaWriter::doOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
			     unsigned long simstep, std::list<ChemicalPotential>* /*lmu*/, map<unsigned, CavityEnsemble>* /*mcav*/ )
{
	domain->calculateGamma(particleContainer,domainDecomp);
	if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		_gammaStream << simstep << "\t"; 
		for (unsigned i=0; i<domain->getNumberOfComponents(); i++){
			_gammaStream << domain->getGamma(i)/_writeFrequency << "\t";
		}
		_gammaStream << endl;
		domain->resetGamma();
	}
}

void GammaWriter::finishOutput(ParticleContainer* /*particleContainer*/,
				DomainDecompBase* /*domainDecomp*/, Domain* /*domain*/){
	_gammaStream.close();
}
