
#include "io/GammaWriter.h"

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "utils/xmlfileUnits.h"
#include "utils/FileUtils.h"

#include <ctime>


using namespace std;

void GammaWriter::readXML(XMLfileUnits& xmlconfig) {
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    global_log->info() << "[GammaWriter] Write frequency: " << _writeFrequency << endl;
    xmlconfig.getNodeValue("outputprefix", _outputPrefix);
    global_log->info() << "[GammaWriter] Output prefix: " << _outputPrefix << endl;

	_range.ymax = global_simulation->getDomain()->getGlobalLength(1);
	xmlconfig.getNodeValue("range/ymin", _range.ymin);
	xmlconfig.getNodeValue("range/ymax", _range.ymax);
	global_log->info() << "[GammaWriter] Range: y: " << _range.ymin << " - " << _range.ymax << endl;
}

void GammaWriter::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	if(domainDecomp->getRank() == 0){
		string resultfilename(_outputPrefix + ".gamma");
		_gammaStream.open(resultfilename);
		_gammaStream.precision(6);
		time_t now;
		time(&now);
		_gammaStream << "# mardyn MD simulation starting at " << ctime(&now) << endl;
		_gammaStream << setw(24) << "simstep" << setw(24) << "gamma" << endl;
	}
}

void GammaWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                          unsigned long simstep)
{
	calculateGamma(particleContainer, domainDecomp);
	if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		double globalLength[3];
		for(int d = 0; d < 3; ++d) {
			globalLength[d] = domain->getGlobalLength(d);
		}
		_gammaStream << FORMAT_SCI_MAX_DIGITS << simstep;
		for(unsigned int componentId = 0; componentId < domain->getNumberOfComponents(); ++componentId){
			_gammaStream << FORMAT_SCI_MAX_DIGITS << getGamma(componentId, globalLength)/_writeFrequency;
		}
		_gammaStream << endl;
		resetGamma();
	}
}

void GammaWriter::finish(ParticleContainer * /*particleContainer*/,
						 DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/){
	_gammaStream.close();
}

void GammaWriter::resetGamma() {
	for(unsigned componentId = 0; componentId < _simulation.getEnsemble()->getComponents()->size(); ++componentId) {
		_Gamma[componentId] = 0;
	}
}

double GammaWriter::getGamma(unsigned id, double globalLength[3]){
	return (_Gamma[id]/(2*globalLength[0]*globalLength[2]));
}

void GammaWriter::calculateGamma(ParticleContainer* particleContainer, DomainDecompBase* domainDecom){
	unsigned numComp = _simulation.getEnsemble()->getComponents()->size();
	std::vector<double> localGamma(numComp, 0.0);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{

		std::vector<double> localGamma_thread(numComp, 0.0);

		for (auto tempMol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {

			// Only for particles within range
			if ((tempMol->r(1) - _range.ymax)*(tempMol->r(1) - _range.ymax) <= 0) {
				unsigned cid = tempMol->componentid();
				localGamma_thread[cid] += tempMol->Vi(1) - 0.5 * (tempMol->Vi(0) + tempMol->Vi(2));
			}
		}

		#if defined(_OPENMP)
		#pragma omp critical
		#endif
		{
			for (unsigned i = 0; i < numComp; ++i) {
				localGamma[i] += localGamma_thread[i];
			}
		}
	}

	domainDecom->collCommInit(numComp);
	for (unsigned i=0; i<numComp; i++){
		domainDecom->collCommAppendDouble(localGamma[i]);
	}
	domainDecom->collCommAllreduceSum();
	for (unsigned i=0; i<numComp; i++){
		localGamma[i] = domainDecom->collCommGetDouble();
	}
	domainDecom->collCommFinalize();
	for (unsigned i=0; i<numComp; i++){
		_Gamma[i]+=localGamma[i];
	}
}
