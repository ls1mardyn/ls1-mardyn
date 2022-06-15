
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

void GammaWriter::init(ParticleContainer */* particleContainer */, DomainDecompBase *domainDecomp, Domain *domain) {
    if(domainDecomp->getRank() == 0){
        string resultfilename(_outputPrefix + ".dat");
        _gammaStream.open(resultfilename);
        _gammaStream.precision(6);
        time_t now;
        time(&now);
        _gammaStream << "# mardyn MD simulation starting at " << ctime(&now) << endl;
        _gammaStream << setw(24) << "simstep";
        for(unsigned int componentId = 0; componentId < domain->getNumberOfComponents() + 1; ++componentId) {
            _gammaStream << setw(22) << "gamma[" << componentId << "]";
        }
        _gammaStream << endl;
        _gammaStream.close();
    }
}

void GammaWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                          unsigned long simstep)
{
    calculateGamma(particleContainer, domainDecomp);

    // Rank 0 writes data to file; write every __writeFrequency steps; do not write data directly after (re)start
    if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0) && (simstep > global_simulation->getNumInitTimesteps())) {
        string resultfilename(_outputPrefix + ".dat");
        double globalLength[3];
        for(int d = 0; d < 3; ++d) {
            globalLength[d] = domain->getGlobalLength(d);
        }
        _gammaStream.open(resultfilename, std::ios::app);
        _gammaStream << FORMAT_SCI_MAX_DIGITS << simstep;
        for(unsigned int componentId = 0; componentId < domain->getNumberOfComponents() + 1; ++componentId) {
            _gammaStream << FORMAT_SCI_MAX_DIGITS << getGamma(componentId, globalLength)/_writeFrequency;
        }
        _gammaStream << endl;
        _gammaStream.close();
        resetGamma(domain);
    }
}

void GammaWriter::resetGamma(Domain *domain) {
    for(unsigned componentId = 0; componentId < domain->getNumberOfComponents() + 1; ++componentId) {
        _Gamma[componentId] = 0;
    }
}

double GammaWriter::getGamma(unsigned id, double globalLength[3]) {
    return (_Gamma[id]/(2*globalLength[0]*globalLength[2]));
}

void GammaWriter::calculateGamma(ParticleContainer* particleContainer, DomainDecompBase* domainDecom) {
    const unsigned numComp = _simulation.getEnsemble()->getComponents()->size() + 1;  // 0 stands for all components
    std::vector<double> localGamma(numComp, 0.0);

    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {

        std::vector<double> localGamma_thread(numComp, 0.0);

        for (auto tempMol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {

            // Only for particles within range
            if ((tempMol->r(1) - _range.ymax)*(tempMol->r(1) - _range.ymin) <= 0) {
                const unsigned cid = tempMol->componentid();
                const double gamma = tempMol->Vi(1) - 0.5 * (tempMol->Vi(0) + tempMol->Vi(2));
                localGamma_thread[cid] += gamma;
                localGamma_thread[0]   += gamma;  // 0 is component-independent value
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