
#include "io/GammaWriter.h"

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "utils/xmlfileUnits.h"

#include <ctime>


using namespace std;

void GammaWriter::readXML(XMLfileUnits& xmlconfig) {
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    global_log->info() << "GammaWriter: Write frequency: " << _writeFrequency << endl;
    xmlconfig.getNodeValue("outputprefix", _outputPrefix);
    global_log->info() << "GammaWriter: Output prefix: " << _outputPrefix << endl;
	xmlconfig.getNodeValue("numInterfaces", _numInterfaces);
    global_log->info() << "GammaWriter: Number of interfaces: " << _numInterfaces << endl;
}

void GammaWriter::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	if(domainDecomp->getRank() == 0){
		string resultfilename(_outputPrefix + ".dat");
		_gammaStream.open(resultfilename);
		_gammaStream.precision(6);
		const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		tm unused{};
		_gammaStream << "# mardyn MD simulation starting at " << std::put_time(localtime_r(&now, &unused), "%c") << endl;
		_gammaStream << "#\tgamma" << endl;
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
		_gammaStream << simstep;
		for(unsigned int componentId = 0; componentId < domain->getNumberOfComponents(); ++componentId){
			_gammaStream << "\t" << getGamma(componentId, globalLength)/_writeFrequency;
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
	// Depending on the number of interfaces
	return (_Gamma[id]/(globalLength[0]*globalLength[2]*_numInterfaces));
}

void GammaWriter::calculateGamma(ParticleContainer* particleContainer, DomainDecompBase* domainDecom){
	unsigned numComp = _simulation.getEnsemble()->getComponents()->size();
	std::vector<double> _localGamma(numComp, 0.0);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{

		std::vector<double> localGamma_thread(numComp, 0.0);

		for (auto tempMol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {

			if (particleContainer->isInBoundingBox(tempMol->r_arr().data())) {
				unsigned cid = tempMol->componentid();
				localGamma_thread[cid] += tempMol->Vi(1) - 0.5 * (tempMol->Vi(0) + tempMol->Vi(2));
			}
		}

		#if defined(_OPENMP)
		#pragma omp critical
		#endif
		{
			for (unsigned i = 0; i < numComp; ++i) {
				_localGamma[i] += localGamma_thread[i];
			}
		}
	}

	domainDecom->collCommInit(numComp);
	for (unsigned i=0; i<numComp; i++){
		domainDecom->collCommAppendDouble(_localGamma[i]);
	}
	domainDecom->collCommAllreduceSum();
	for (unsigned i=0; i<numComp; i++){
		_localGamma[i] = domainDecom->collCommGetDouble();
	}
	domainDecom->collCommFinalize();
	for (unsigned i=0; i<numComp; i++){
		_Gamma[i]+=_localGamma[i];
	}
}
