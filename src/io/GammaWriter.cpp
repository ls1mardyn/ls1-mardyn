
#include "io/GammaWriter.h"

#include <algorithm>
#include <vector>

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"
#include "Simulation.h"
#include "utils/xmlfileUnits.h"
#include "utils/FileUtils.h"


void GammaWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "[GammaWriter] Write frequency: " << _writeFrequency << std::endl;
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "[GammaWriter] Output prefix: " << _outputPrefix << std::endl;
	xmlconfig.getNodeValue("numInterfaces", _numInterfaces);
	Log::global_log->info() << "[GammaWriter] Number of interfaces: " << _numInterfaces << std::endl;

	_range.ymax = global_simulation->getDomain()->getGlobalLength(1);
	xmlconfig.getNodeValue("range/ymin", _range.ymin);
	xmlconfig.getNodeValue("range/ymax", _range.ymax);
	Log::global_log->info() << "[GammaWriter] Range: y: " << _range.ymin << " - " << _range.ymax << std::endl;
}

void GammaWriter::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {

	for (int d = 0; d < 3; ++d) {
		_globalLength[d] = domain->getGlobalLength(d);
	}
	_numComp = domain->getNumberOfComponents() + 1;  // 0 stands for all components

	_gamma.resize(_numComp);

	// Rank 0 writes data to file
	if (domainDecomp->getRank() == 0) {
		std::string resultfilename(_outputPrefix + ".dat");
		_gammaStream.open(resultfilename);
		_gammaStream.precision(6);
		_gammaStream << std::setw(24) << "simstep";
		for (unsigned int componentId = 0; componentId < _numComp; ++componentId) {
			_gammaStream << std::setw(22) << "gamma[" << componentId << "]";
		}
		_gammaStream << std::endl;
		_gammaStream.close();
	}
}

void GammaWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
						  unsigned long simstep)
{
	calculateGamma(particleContainer, domainDecomp);

	// Write every _writeFrequency steps; do not write data directly after (re)start
	if ((simstep % _writeFrequency == 0) && (simstep > global_simulation->getNumInitTimesteps())) {
		// Rank 0 writes data to file
		if (domainDecomp->getRank() == 0) {
			std::string resultfilename(_outputPrefix + ".dat");

			_gammaStream.open(resultfilename, std::ios::app);
			_gammaStream << FORMAT_SCI_MAX_DIGITS << simstep;
			for (unsigned int componentId = 0; componentId < _numComp; ++componentId) {
				_gammaStream << FORMAT_SCI_MAX_DIGITS << getGamma(componentId)/_writeFrequency;
			}
			_gammaStream << std::endl;
			_gammaStream.close();
		}
		resetGamma(domain);
	}
}

inline void GammaWriter::resetGamma(Domain *domain) {
	for (unsigned componentId = 0; componentId < _numComp; ++componentId) {
		_gamma[componentId] = 0.0;
	}
}

inline double GammaWriter::getGamma(unsigned id) {
	return (_gamma[id]/(_globalLength[0]*_globalLength[2]*_numInterfaces));
}

void GammaWriter::calculateGamma(ParticleContainer* particleContainer, DomainDecompBase* domainDecom) {
	std::vector<double> localGamma(_numComp, 0.0);

	{
	#if defined(_OPENMP)
	// Taken from https://stackoverflow.com/a/43169193 and modified
	#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
								std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
								initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
	#pragma omp parallel reduction(vec_double_plus: localGamma)
	#endif
		for (auto tempMol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {
			// Only for particles within range
			if ((tempMol->r(1) - _range.ymax)*(tempMol->r(1) - _range.ymin) <= 0) {
				const unsigned cid = tempMol->componentid() + 1;
				const double gamma = tempMol->Vi(1) - 0.5 * (tempMol->Vi(0) + tempMol->Vi(2));
				localGamma[cid] += gamma;
				localGamma[0]   += gamma;  // 0 is component-independent value
			}
		}
	}

	domainDecom->collCommInit(_numComp);
	for (unsigned int i=0; i<_numComp; i++) {
		domainDecom->collCommAppendDouble(localGamma[i]);
	}
	domainDecom->collCommAllreduceSum();
	for (unsigned int i=0; i<_numComp; i++) {
		localGamma[i] = domainDecom->collCommGetDouble();
	}
	domainDecom->collCommFinalize();
	for (unsigned int i=0; i<_numComp; i++) {
		_gamma[i] += localGamma[i];
	}
}
