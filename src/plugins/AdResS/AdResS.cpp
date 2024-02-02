//
// Created by alex on 27.04.23.
//

#include "AdResS.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "Domain.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "plugins/AdResS/parallel/AdResSKDDecomposition.h"
#include "utils/mardyn_assert.h"
#include "Interpolation.h"

#include <cmath>

using namespace std;
using Log::global_log;

AdResS::AdResS() : _fpRegions(), _particleContainer(nullptr),
                   _components(nullptr), _comp_to_res(), _domain(nullptr) {};

AdResS::~AdResS() = default;

void AdResS::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
    global_log->debug() << "[AdResS] Enabled " << std::endl;
    for(const auto& region : _fpRegions) {
        global_log->debug() << "[AdResS] FPRegion Box from ["
        << region._low[0] << "," << region._low[1] << "," << region._low[2] << "] to ["
        << region._high[0] << "," << region._high[1] << "," << region._high[2] << "] with hybrid dims ["
        << region._hybridDims[0] << "," << region._hybridDims[1] << "," << region._hybridDims[2] << "]" << std::endl;
    }

    _components = _simulation.getEnsemble()->getComponents();
    _domain = domain;
    _particleContainer = particleContainer;
    _comp_to_res.resize(_components->size());

    for(Component& comp : *_components) {
        unsigned int id = comp.ID();
        if(id % Resolution::ResolutionCount == Resolution::FullParticle) {
            if(comp.getName().rfind("FP_", 0) == 0) _comp_to_res[id] = FullParticle;
            else {
                global_log->fatal() << "[AdResS] Component with id=" << id+1 << " does not have FP_ prefix in name."
                << "Uncertain if this component should be Full Particle representation." << std::endl;
                exit(669);
            }
        }
        if(id % Resolution::ResolutionCount == Resolution::Hybrid) {
            if(comp.getName().rfind("H_", 0) == 0) _comp_to_res[id] = Hybrid;
            else {
                global_log->fatal() << "[AdResS] Component with id=" << id+1 << " does not have H_ prefix in name."
                                    << "Uncertain if this component should be Hybrid representation." << std::endl;
                exit(669);
            }
        }
        if(id % Resolution::ResolutionCount == Resolution::CoarseGrain) {
            if(comp.getName().rfind("CG_", 0) == 0) _comp_to_res[id] = CoarseGrain;
            else {
                global_log->fatal() << "[AdResS] Component with id=" << id+1 << " does not have CG_ prefix in name."
                                    << "Uncertain if this component should be Coarse Grain representation." << std::endl;
                exit(669);
            }
        }
    }
    if(_enableThermodynamicForce) {
        _densityProfiler.init(_samplingStepSize, domain, _rho0, _smoothingFactor);
        _thermodynamicForceSampleCounter = 0;
        _thermodynamicForceHist.n = 0;
    }
    if(_enableThermodynamicForce && _createThermodynamicForce) {
        _densityProfiler.sampleDensities(particleContainer, domainDecomp, domain);
        _targetDensity = std::vector<double>(_densityProfiler.getDensity(0));
        _thermodynamicForce.n = _targetDensity.size();
        _thermodynamicForce.begin = 0.0;
        _thermodynamicForce.step_width.resize(_thermodynamicForce.n-1, _samplingStepSize);
        _thermodynamicForce.gradients.resize(_thermodynamicForce.n, 0.0);
        _thermodynamicForce.function_values.resize(_thermodynamicForce.n, 0.0);

        _thermodynamicForceHist.n = _targetDensity.size();
        _thermodynamicForceHist.begin = 0.0;
        _thermodynamicForceHist.step_width.resize(_thermodynamicForceHist.n-1, _samplingStepSize);
        _thermodynamicForceHist.gradients.resize(_thermodynamicForceHist.n, 0.0);
        _thermodynamicForceHist.function_values.resize(_thermodynamicForceHist.n, 0.0);

        if(_logDensities) writeDensities("F_TH_TargetDensity.txt", _targetDensity);
    }

    _tracerProcessor = new TracerCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), _simulation.getParticlePairsHandler(), _comp_to_res);
    _simulation.setCellProcessor(_tracerProcessor);
}

void AdResS::readXML(XMLfileUnits &xmlconfig) {
    XMLfile::Query query = xmlconfig.query("enableFTH");
    _enableThermodynamicForce = false;
    unsigned long count = query.card();
    if (count > 1) {
        global_log->fatal() << "[AdResS] Only specify one enableFTH block in config file!" << std::endl;
        Simulation::exit(668);
    }
    //F_TH enabled
    if(count == 1) {
        _enableThermodynamicForce = true;
        _forceMax = xmlconfig.getNodeValue_double("enableFTH/fmax", 100.0);
        _forceMax = std::abs(_forceMax);
        _logFTH = xmlconfig.getNodeValue_bool("enableFTH/logFTH", false);
        _logDensities = xmlconfig.getNodeValue_bool("enableFTH/logDensity", false);
        _rho0 = xmlconfig.getNodeValue_double("enableFTH/rho", 0.01);
        _smoothingFactor = xmlconfig.getNodeValue_double("enableFTH/smoothingFactor", 10);

        query = xmlconfig.query("enableFTH/createFTH");
        count = query.card();
        if (count > 1) {
            global_log->fatal() << "[AdResS] Only specify one createFTH block in config file!" << std::endl;
            Simulation::exit(668);
        }

        if(count == 1) { // sample FTH function
            _createThermodynamicForce = true;
            _thermodynamicForceSampleGap = xmlconfig.getNodeValue_int("enableFTH/createFTH/sampleGap", 200);
            _convergenceThreshold = xmlconfig.getNodeValue_double("enableFTH/createFTH/threshold", 0.02);
            _convergenceFactor = xmlconfig.getNodeValue_double("enableFTH/createFTH/convFactor", 0.2);
            _samplingStepSize = xmlconfig.getNodeValue_double("enableFTH/createFTH/sampleBinSize", 0.2);

        }
        else { // use existing FTH function
            query = xmlconfig.query("enableFTH/forceFunction");
            count = query.card();
            if (count != 1) {
                global_log->fatal() << "[AdResS] Must specify one forceFunction block in config file!" << std::endl;
                Simulation::exit(668);
            }
            _logFTH = false;
            _thermodynamicForce.begin = xmlconfig.getNodeValue_double("enableFTH/forceFunction/startX");
            query = xmlconfig.query("enableFTH/forceFunction/samplePoint");
            count = query.card();
            if(count < 5) {
                global_log->fatal() << "[AdResS] Force function must have at least 5 sample points!" << std::endl;
                Simulation::exit(668);
            }

            _thermodynamicForce.n = count;
            _thermodynamicForce.gradients.resize(count, 0.0);
            _thermodynamicForce.function_values.resize(count, 0.0);
            _thermodynamicForce.step_width.resize(count-1, 0.0);
            XMLfile::Query::const_iterator sampleIter;
            std::string oldpath = xmlconfig.getcurrentnodepath();
            for(sampleIter = query.begin(); sampleIter; sampleIter++) {
                xmlconfig.changecurrentnode(sampleIter);
                unsigned long id = 0;
                xmlconfig.getNodeValue("@id", id);
                xmlconfig.getNodeValue("grad", _thermodynamicForce.gradients[id - 1]);
                xmlconfig.getNodeValue("func", _thermodynamicForce.function_values[id - 1]);
                if(id < count) xmlconfig.getNodeValue("step", _thermodynamicForce.step_width[id - 1]);
            }
            xmlconfig.changecurrentnode(oldpath);
            _samplingStepSize = _thermodynamicForce.step_width[0];
            _thermodynamicForceSampleGap = 200;
            _createThermodynamicForce = false;
        }
    }

    long numRegions = 0;
    query = xmlconfig.query("fpregions/region");
    numRegions = query.card();
    if (numRegions == 0) {
        global_log->fatal() << "No AdResS regions specified: Falling back to dynamic region selection. ERROR: not implemented yet!" << std::endl;
        Simulation::exit(668);
    }
    if (numRegions != 1) {
        global_log->fatal() << "AdResS currently only supports a single FP Region. For more general support select a previous version." << std::endl;
        Simulation::exit(-1);
    }
    _fpRegions.resize(numRegions);

    XMLfile::Query::const_iterator regionIter;
    std::string oldpath = xmlconfig.getcurrentnodepath();
    for(regionIter = query.begin(); regionIter; regionIter++) {
        xmlconfig.changecurrentnode(regionIter);
        unsigned int id = 0;
        xmlconfig.getNodeValue("@id", id);
        _fpRegions[id - 1].readXML(xmlconfig);
    }
    xmlconfig.changecurrentnode(oldpath);

    _domain = _simulation.getDomain();
    _comp_to_res.resize(_simulation.getEnsemble()->getComponents()->size(), FullParticle);
    if(auto decomp = dynamic_cast<AdResSKDDecomposition*>(&_simulation.domainDecomposition())) {
        decomp->setAdResSPlugin(this);
    }

    for(const auto& region : _fpRegions) {
        global_log->info() << "[AdResS] FPRegion Box from ["
                            << region._low[0] << "," << region._low[1] << "," << region._low[2] << "] to ["
                            << region._high[0] << "," << region._high[1] << "," << region._high[2] << "] with hybrid dim "
                            << region._hybridDims[0] << "," << region._hybridDims[1] << "," << region._hybridDims[2] << "]" << std::endl;
    }
}

void AdResS::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                     unsigned long simstep) {
    if(_thermodynamicForceSampleCounter != 0) return;

    std::stringstream stream;
    stream << "./F_TH_InterpolationFunction_" << simstep << ".xml";
    if(_logFTH) writeFunctionToXML(stream.str(), _thermodynamicForce);

    stream.clear();
    stream = std::stringstream {};
    stream << "./F_TH_Density_" << simstep << ".txt";
    _densityProfiler.sampleDensities(particleContainer, domainDecomp, domain);
    std::vector<double> densities{_densityProfiler.getDensity(0)};
    if(_logDensities) writeDensities(stream.str(), densities);

    stream.clear();
    stream = std::stringstream {};
    stream << "./F_TH_Density_Smooth_" << simstep << ".txta";
    densities = _densityProfiler.getDensitySmoothed(0);
    if(_logDensities) writeDensities(stream.str(), densities);

    _densityProfiler.computeDensities(particleContainer, domainDecomp, domain);
    Interpolation::Function histDensity = _densityProfiler.getHistDensity(0);
    Interpolation::resampleFunction(0, domain->getGlobalLength(0), _samplingStepSize, histDensity);
    stream.clear();
    stream = std::stringstream {};
    stream << "./F_TH_HistDensity_" << simstep << ".xmla";
    writeFunctionToXML(stream.str(), histDensity);

    stream.clear();
    stream = std::stringstream {};
    stream << "./F_TH_HistInterpolationFunction_" << simstep << ".xmlb";
    writeFunctionToXML(stream.str(), _thermodynamicForceHist);
}

void AdResS::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {

}

std::string AdResS::getPluginName() {
    return {"AdResS"};
}

void AdResS::beforeForces(ParticleContainer *container, DomainDecompBase *, unsigned long) {
    // check for all particles if it is in a certain region to set component correctly
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    for(auto itM = container->iterator(ParticleIterator::ALL_CELLS); itM.isValid(); ++itM) {
        bool stop = false;
        //check if is in any full particle region
        for(auto& reg : _fpRegions) {
            if(reg.isInnerPointDomain(_domain, FullParticle, itM->r_arr())) {
                checkMoleculeLOD(*itM, FullParticle);
                stop = true;
                break;
            }
        }
        if(stop) continue;

        //check if is in any hybrid region
        for(auto& reg : _fpRegions) {
            if(reg.isInnerPointDomain(_domain, Hybrid, itM->r_arr())) {
                checkMoleculeLOD(*itM, Hybrid);
                stop = true;
                break;
            }
        }
        if(stop) continue;

        // molecule is in no Hybrid or FP region
        checkMoleculeLOD(*itM, CoarseGrain);
    }

    // handle thermodynamic force
    if(!_enableThermodynamicForce) return;
    _thermodynamicForceSampleCounter++;

    if(_thermodynamicForceSampleCounter % _thermodynamicForceSampleGap != 0) return;
    _thermodynamicForceSampleCounter = 0;

    if(!_createThermodynamicForce) return;
    if(checkF_TH_Convergence()) {
        global_log->info() << "[AdResS] F_TH has converged." << std::endl;
        writeFunctionToXML("./F_TH_InterpolationFunction_Final.xml", _thermodynamicForce);
        Simulation::exit(0);
    }
    else computeF_TH();
}

void AdResS::siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
    if(_enableThermodynamicForce) {
        applyF_TH();
    }
}

void AdResS::computeF_TH() {
    _densityProfiler.sampleDensities(_particleContainer, &_simulation.domainDecomposition(), _simulation.getDomain());
    std::vector<double> d{_densityProfiler.getDensitySmoothed(0)};
    std::vector<double> d_prime;
    Interpolation::computeGradient(d, d_prime);
    Interpolation::Function d_prime_fun;
    std::vector<double> steps;
    steps.resize(d_prime.size()-1, _samplingStepSize);
    Interpolation::computeHermite(0.0, d_prime, steps, d_prime.size(), d_prime_fun);

    _densityProfiler.computeDensities(_particleContainer, &_simulation.domainDecomposition(), _simulation.getDomain());
    Interpolation::Function d_hist = _densityProfiler.getHistDensity(0);
    std::vector<double> fVals;
    fVals.resize(_thermodynamicForceHist.n);
    double pos = 0.0;
    for(int i = 0; i < _thermodynamicForceHist.n; i++) {
        fVals[i] = Interpolation::computeHermiteAt(pos, d_hist);
        pos += _thermodynamicForceHist.step_width[i];
    }
    std::vector<double> steps_hist;
    steps_hist.resize(_thermodynamicForceHist.n-1, _samplingStepSize);
    Interpolation::computeHermite(0.0, fVals, steps_hist, _thermodynamicForceHist.n, d_hist);
    Interpolation::Function d_prime_hist;
    Interpolation::computeHermite(0.0, d_hist.gradients, d_hist.step_width, d_hist.n, d_prime_hist);

    for(int i = 0; i < d_prime_fun.n; i++) {
        _thermodynamicForce.function_values[i] -= _convergenceFactor * d_prime_fun.function_values[i];
        _thermodynamicForce.gradients[i] -= _convergenceFactor * d_prime_fun.gradients[i];

        _thermodynamicForceHist.function_values[i] -= _convergenceFactor * d_prime_hist.function_values[i];
        _thermodynamicForceHist.gradients[i] -= _convergenceFactor * d_prime_hist.gradients[i];
    }

    // TODO FIXME!!
    auto& region = _fpRegions[0];
    double x_pos = _thermodynamicForce.begin;
    double low = region._low[0];
    double high = region._high[0];
    for(unsigned long i = 0; i < _thermodynamicForce.n; i++) {
        if(x_pos >= low && x_pos <= high){
            _thermodynamicForce.function_values[i] = 0;
            _thermodynamicForce.gradients[i] = 0;

            _thermodynamicForceHist.function_values[i] = 0;
            _thermodynamicForceHist.gradients[i] = 0;
        }
        x_pos += _thermodynamicForce.step_width[i];
    }

    _lastGradient = std::move(d_prime_fun);
}

bool AdResS::checkF_TH_Convergence() {
    if(_lastGradient.function_values.empty()) return false;
    auto it0 = std::max_element(_lastGradient.function_values.begin(), _lastGradient.function_values.end());
    auto it1 = std::min_element(_lastGradient.function_values.begin(), _lastGradient.function_values.end());
    double max_grad = std::max(std::abs(*it0), std::abs(*it1)) * _convergenceFactor;
    global_log->info() << "[AdResS] F_TH conv delta: " << max_grad << std::endl;
    return max_grad <= _convergenceThreshold;
}

void AdResS::applyF_TH() {
    std::array<double, 3> low = {2*_samplingStepSize,
                                 0,
                                 0};
    std::array<double, 3> high= {_simulation.getDomain()->getGlobalLength(0) - 2*_samplingStepSize,
                                 _simulation.getDomain()->getGlobalLength(1),
                                 _simulation.getDomain()->getGlobalLength(2)};
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    for (auto itM = _particleContainer->regionIterator(std::data(low), std::data(high), ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
        if(_comp_to_res[itM->componentid()] == FullParticle) continue;

        double x = itM->r(0);
        double F = computeHermiteAt(x, _thermodynamicForce);
        std::array<double, 3> force = {F, 0.0, 0.0};
        itM->Fadd(std::data(force));
    }
}

void AdResS::checkMoleculeLOD(Molecule &molecule, Resolution targetRes) {
    auto id = molecule.component()->ID();
    //molecule has already correct component
    if(_comp_to_res[id] == targetRes) return;

    int offset = static_cast<int>(targetRes) - static_cast<int>(_comp_to_res[id]);
    molecule.setComponent(&_components->at(id + offset));
}

void AdResS::writeFunctionToXML(const string &filename, Interpolation::Function &fun) {
#ifdef ENABLE_MPI
    if (_simulation.domainDecomposition().getRank() != 0) return;
#endif
    try {
        std::ofstream file {filename};
        file << "<forceFunction>\n";
        file << "    <startX>" << fun.begin << "</startX>\n";
        for(unsigned long i = 0; i < fun.n; i++) {
            file << "    <samplePoint id=\"" << i + 1 <<"\">\n";
            file << "    " << "    <grad>" << fun.gradients[i] << "</grad>\n";
            file << "    " << "    <func>" << fun.function_values[i] << "</func>\n";
            if(i < fun.n - 1) file << "    " << "    <step>" << fun.step_width[i] << "</step>\n";
            file << "    </samplePoint>\n";
        }
        file << "</forceFunction>\n";
        file.flush();
        file.close();
    } catch (std::ifstream::failure& e) {
        global_log->error() << "[AdResS] Failed to write Interpolation function.\n" << e.what() << std::endl;
        _simulation.exit(-1);
    }
}

void AdResS::writeDensities(const string &filename, vector<double> &densities, const string &separator) {
#ifdef ENABLE_MPI
    if (_simulation.domainDecomposition().getRank() != 0) return;
#endif
    try {
        std::ofstream file {filename};
        for(unsigned long i = 0; i < densities.size()-1; i++) {
            file << densities[i] << separator;
        }
        file << densities[densities.size()-1] << std::endl;
    } catch (std::ifstream::failure& e) {
        global_log->error() << "[AdResS] Failed to write densities to file.\n" << e.what() << std::endl;
        _simulation.exit(-1);
    }
}

void AdResS::afterForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
    std::array<double, 3> low = {2*_samplingStepSize,
                                 0,
                                 0};
    std::array<double, 3> high= {_simulation.getDomain()->getGlobalLength(0) - 2*_samplingStepSize,
                                 _simulation.getDomain()->getGlobalLength(1),
                                 _simulation.getDomain()->getGlobalLength(2)};
    double cutoff = _simulation.getcutoffRadius();
    auto& region = _fpRegions[0];
    low[0] = region._lowHybrid[0] - cutoff;
    high[0] = region._low[0] + cutoff;
#if defined(_OPENMP)
#pragma omp parallel
#endif
    for (auto itM = _particleContainer->regionIterator(std::data(low), std::data(high), ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
        std::array<double, 3> f = itM->F_arr();
        for(short d = 0; d < 3; d++) {
            f[d] = std::copysign(std::min(std::abs(f[d]), _forceMax), f[d]);
        }
        itM->setF(std::data(f));
    }

    low[0] = region._high[0] - cutoff;
    high[0] = region._highHybrid[0] + cutoff;
#if defined(_OPENMP)
#pragma omp parallel
#endif
    for (auto itM = _particleContainer->regionIterator(std::data(low), std::data(high), ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
        std::array<double, 3> f = itM->F_arr();
        for(short d = 0; d < 3; d++) {
            f[d] = std::copysign(std::min(std::abs(f[d]), _forceMax), f[d]);
        }
        itM->setF(std::data(f));
    }
}
