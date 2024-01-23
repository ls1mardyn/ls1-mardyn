//
// Created by alex on 27.04.23.
//

#include "AdResS.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "Domain.h"
#include "plugins/AdResS/util/AdResSRegionTraversal.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "plugins/AdResS/parallel/AdResSKDDecomposition.h"
#include "utils/mardyn_assert.h"
#include "Interpolation.h"

#include <cmath>

using namespace std;
using Log::global_log;

double (*AdResS::weight)(std::array<double, 3> r, FPRegion& region) = nullptr;

AdResS::AdResS() : _mesoVals(), _forceAdapter(nullptr), _fpRegions(), _particleContainer(nullptr),
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

    if(_enableThermodynamicForce && _createThermodynamicForce) {
        _densityProfiler.init(_samplingStepSize, domain, _rho0);
        _densityProfiler.sampleDensities(particleContainer, domainDecomp, domain);
        _targetDensity = std::vector<double>(_densityProfiler.getDensity(0));
        _thermodynamicForceSampleCounter = 0;
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
}

void AdResS::readXML(XMLfileUnits &xmlconfig) {
    static const std::array<std::string, 5> weight_impls {"euclid", "manhattan", "component", "near", "flat"};
    static const std::array<double (*)(std::array<double, 3> r, FPRegion &region), 5> impls
        {weightEuclid, weightManhattan, weightComponent, weightNearest, weightFlat};
    std::string impl = weight_impls[0];
    xmlconfig.getNodeValue("weightImpl", impl);
    unsigned long index = std::distance(weight_impls.begin(), std::find(weight_impls.begin(), weight_impls.end(), impl));
    AdResS::weight = impls[index >= 5 ? 0 : index];
    global_log->info() << "[AdResS] Using weight implementation " << weight_impls[index >= 5 ? 0 : index] << std::endl;

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

        query = xmlconfig.query("enableFTH/createFTH");
        count = query.card();
        if (count > 1) {
            global_log->fatal() << "[AdResS] Only specify one createFTH block in config file!" << std::endl;
            Simulation::exit(668);
        }

        if(count == 1) { // sample FTH function
            _createThermodynamicForce = true;
            _thermodynamicForceSampleGap = xmlconfig.getNodeValue_int("enableFTH/createFTH/sampleGap", 100) - 1;
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
    // todo add check that no region overlap even in hybrid considering periodic bounds

    _forceAdapter = new AdResSForceAdapter(*this);
    _simulation.setParticlePairsHandler(_forceAdapter);
    _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), _forceAdapter));
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
    if(_enableThermodynamicForce) {
        if(_createThermodynamicForce && _thermodynamicForceSampleCounter >= _thermodynamicForceSampleGap) {
            if(checkF_TH_Convergence()) {
                global_log->info() << "[AdResS] F_TH has converged." << std::endl;
                writeFunctionToXML("./F_TH_InterpolationFunction_Final.xml", _thermodynamicForce);
                Simulation::exit(0);
            }
            else computeF_TH();
            _thermodynamicForceSampleCounter = -1;
        }
        _thermodynamicForceSampleCounter++;
    }
}

void AdResS::siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
    //computeForce(true);
    //computeForce(false);
    _mesoVals.setInDomain(_domain);
    _mesoVals.clear();

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
                                 _simulation.getDomain()->getGlobalLength(1),
                                 _simulation.getDomain()->getGlobalLength(2)};
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

    // TODO FIXME!!
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

void AdResS::checkMoleculeLOD(Molecule &molecule, Resolution targetRes) {
    auto id = molecule.component()->ID();
    //molecule has already correct component
    if(_comp_to_res[id] == targetRes) return;

    int offset = static_cast<int>(targetRes) - static_cast<int>(_comp_to_res[id]);
    molecule.setComponent(&_components->at(id + offset));
}

void AdResS::computeForce(bool invert) {
    double cutoff = _simulation.getcutoffRadius();
    double cutoff2 = cutoff * cutoff;
    double LJCutoff2 = _simulation.getLJCutoff();
    LJCutoff2 *= LJCutoff2;
    std::array<double,3> globLen{0};
    std::array<double, 3> pc_low{0};
    std::array<double, 3> pc_high{0};
    for(int d = 0; d < 3; d++) {
        globLen[d] = _domain->getGlobalLength(d);
        pc_low[d] = _particleContainer->getBoundingBoxMin(d);
        pc_high[d] = _particleContainer->getBoundingBoxMax(d);
    }

    //check all regions
    for(auto& region : _fpRegions) {
        //check for local node if region is even in particle container
        if(!region.isRegionInBox(pc_low, pc_high)) continue;

        // first check if the region crosses any domain bounds
        std::array<bool, 3> inDim {false};
        for(int d = 0; d < 3; d++) {
            inDim[d] = region._lowHybrid[d] >= 0 && region._highHybrid[d] <= globLen[d];
        }
        // find how much it crosses bounds
        std::array<std::array<double,2>,3> deltaLowHighDim{};
        for(int d = 0; d < 3; d++) {
            deltaLowHighDim[d][0] = std::max(-(region._lowHybrid[d]), 0.);
            deltaLowHighDim[d][1] = std::max((region._highHybrid[d]) - globLen[d], 0.);
        }
        // check all regions that wrap around due to periodic bounds
        // we only need to check for each dim twice if it actually wraps around in that dimension
        // in the arrays we create the indices of the boxes depending on the viewed case
        // for x y or z: if they are 0 then we do not check a wrap around in that dimension
        for(int x = 0; x <= !inDim[0]; x++) {
            for(int y = 0; y <= !inDim[1]; y++) {
                for(int z = 0; z <= !inDim[2]; z++) {
                    std::array<double, 3> checkLow {
                            (1 - x) * std::max(region._lowHybrid[0], 0.) + x * ((deltaLowHighDim[0][0] != 0) * (globLen[0] - deltaLowHighDim[0][0]) + (deltaLowHighDim[0][1] != 0) * 0),
                            (1 - y) * std::max(region._lowHybrid[1], 0.) + y * ((deltaLowHighDim[1][0] != 0) * (globLen[1] - deltaLowHighDim[1][0]) + (deltaLowHighDim[1][1] != 0) * 0),
                            (1 - z) * std::max(region._lowHybrid[2], 0.) + z * ((deltaLowHighDim[2][0] != 0) * (globLen[2] - deltaLowHighDim[2][0]) + (deltaLowHighDim[2][1] != 0) * 0) };

                    std::array<double, 3> checkHigh {
                            (1 - x) * std::min(region._highHybrid[0], globLen[0]) + x * ((deltaLowHighDim[0][0] != 0) * (globLen[0]) + (deltaLowHighDim[0][1] != 0) * deltaLowHighDim[0][1]),
                            (1 - y) * std::min(region._highHybrid[1], globLen[1]) + y * ((deltaLowHighDim[1][0] != 0) * (globLen[1]) + (deltaLowHighDim[1][1] != 0) * deltaLowHighDim[1][1]),
                            (1 - z) * std::min(region._highHybrid[2], globLen[2]) + z * ((deltaLowHighDim[2][0] != 0) * (globLen[2]) + (deltaLowHighDim[2][1] != 0) * deltaLowHighDim[2][1]) };
                    // checkLow and checkHigh create a box within bounds, that was potentially wrapped around due to periodic bounds
                    // now create bigger box surrounding this box to find molecules that have interacted with our hybrid molecules
                    // only molecules within cutoff around the region have interacted with hybrid molecules
                    for(int d = 0; d < 3; d++) {
                        checkLow[d] -= cutoff;
                        checkHigh[d] += cutoff;
                    }

                    // now have created a box in which forces need to be calculated
                    // this we will multi-thread: for that we implement a simplified C08-Traversal
                    _forceAdapter->init(_domain); // clear thread data
                    AdResSRegionTraversal traversal{ checkLow, checkHigh, _particleContainer, _comp_to_res};
                    traversal.traverse(*_forceAdapter, region, invert);
                    _forceAdapter->finish(); // gather thread data
                }
            }
        }
    }
}

double AdResS::weightEuclid(std::array<double, 3> r, FPRegion& region) {
    // if point is in the FP region -> weight is 1
    if(FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
    // point is in hybrid region -> weight is between 0 and 1
    else if(FPRegion::isInnerPoint(r, region._lowHybrid, region._highHybrid)) {
        std::array<double, 3> intersect_inner = region.computeIntersection(r, FPRegion::Intersection::H_FP);
        std::array<double, 3> intersect_outer = region.computeIntersection(r, FPRegion::Intersection::CG_H);
        double hyb_axis_length = sqrt(std::pow(intersect_outer[0]-intersect_inner[0], 2) +
                                      std::pow(intersect_outer[1]-intersect_inner[1], 2) +
                                      std::pow(intersect_outer[2]-intersect_inner[2], 2));
        double dist = sqrt(std::pow(r[0]-intersect_inner[0], 2) +
                           std::pow(r[1]-intersect_inner[1], 2) +
                           std::pow(r[2]-intersect_inner[2], 2));

        return std::pow(std::cos(M_PI/(2*hyb_axis_length) * dist), 2);
    }
    // point is in the CG region -> weight is 0
    else return 0.;
}

double AdResS::weightManhattan(std::array<double, 3> r, FPRegion &region) {
    // if point is in the FP region -> weight is 1
    if(FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
    // point is in hybrid region -> weight is between 0 and 1
    else if(FPRegion::isInnerPoint(r, region._lowHybrid, region._highHybrid)) {
        std::array<double, 3> intersect_inner = region.computeIntersection(r, FPRegion::Intersection::H_FP);
        std::array<double, 3> intersect_outer = region.computeIntersection(r, FPRegion::Intersection::CG_H);
        double hyb_axis_length = intersect_outer[0]-intersect_inner[0]  +
                                 intersect_outer[1]-intersect_inner[1]  +
                                 intersect_outer[2]-intersect_inner[2];
        double dist = r[0]-intersect_inner[0] +
                      r[1]-intersect_inner[1] +
                      r[2]-intersect_inner[2];

        return std::pow(std::cos(M_PI/(2*hyb_axis_length) * dist), 2);
    }
    // point is in the CG region -> weight is 0
    else return 0.;
}

double AdResS::weightComponent(std::array<double, 3> r, FPRegion &region) {
    // if point is in the FP region -> weight is 1
    if(FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
    // point is in hybrid region -> weight is between 0 and 1
    else if(FPRegion::isInnerPoint(r, region._lowHybrid, region._highHybrid)) {
        std::array<bool, 3> contribute_dir = {r[0] <= region._low[0] || r[0] >= region._high[0],
                                              r[1] <= region._low[1] || r[1] >= region._high[1],
                                              r[2] <= region._low[2] || r[2] >= region._high[2]};
        std::array<double,3> dist_dir{std::max(std::max(region._low[0] - r[0], 0.), std::max(r[0] - region._high[0], 0.)),
                                      std::max(std::max(region._low[1] - r[1], 0.), std::max(r[1] - region._high[1], 0.)),
                                      std::max(std::max(region._low[2] - r[2], 0.), std::max(r[2] - region._high[2], 0.))};
        double w = 1.;
        for(int d = 0; d < 3; d++) {
            if(region._hybridDims[d] != 0)
                w *= contribute_dir[d] * std::cos(M_PI/(2*region._hybridDims[d]) * dist_dir[d]) + (1 - contribute_dir[d]);
        }
        return w;
    }
    // point is in the CG region -> weight is 0
    else return 0.;
}

double AdResS::weightNearest(std::array<double, 3> r, FPRegion &region) {
    // if point is in the FP region -> weight is 1
    if(FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
    // point is in hybrid region -> weight is between 0 and 1
    else if(FPRegion::isInnerPoint(r, region._lowHybrid, region._highHybrid)) {
        //example: dist_dir = {dx, 0, 0} if point is either left or right of region
        //         dist_dir = {dx, dy, 0} if point is in edged that go front to rear in block
        //         dist_dir = {dx, dy, dz} if point is in one corner
        std::array<double,3> dist_dir{std::max(std::max(region._low[0] - r[0], 0.), std::max(r[0] - region._high[0], 0.)),
                                      std::max(std::max(region._low[1] - r[1], 0.), std::max(r[1] - region._high[1], 0.)),
                                      std::max(std::max(region._low[2] - r[2], 0.), std::max(r[2] - region._high[2], 0.))};
        // => distance to the closest point of r on the surface of the inner region is the l2-norm of dist_dir
        double dist = sqrt(std::pow(dist_dir[0],2)+std::pow(dist_dir[1],2)+std::pow(dist_dir[2],2));
        // compute hybrid size as hybrid width is not equal on all sides
        double hDim = 0;
        for(int d = 0; d < 3; d++) {
            hDim += (dist_dir[d] != 0) * std::pow(region._hybridDims[d], 2);
        }
        hDim = sqrt(hDim);
        if(dist >= hDim) return 0.; // we are outside the rounded region -> treat as CG
        return std::pow(std::cos(M_PI/(2*hDim) * dist), 2);
    }
    // point is in the CG region -> weight is 0
    else return 0.;
}

double AdResS::weightFlat(std::array<double, 3> r, FPRegion &region) {
    if(FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
    return 0.;
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
