//
// Created by alex on 27.04.23.
//

#include "AdResS.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "Domain.h"

#include <cmath>

using namespace std;
using Log::global_log;

AdResS::AdResS() : _mesoVals(), _forceAdapter(_mesoVals), _fpRegions(), _particleContainer(nullptr),
                   _components(nullptr), _comp_to_res(), _domain(nullptr) {};

AdResS::~AdResS() = default;

void AdResS::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
    global_log->debug() << "[AdResS] Enabled " << std::endl;
    for(const auto& region : _fpRegions) {
        global_log->debug() << "[AdResS] FPRegion Box from ["
        << region._low[0] << "," << region._low[1] << "," << region._low[2] << "] to ["
        << region._high[0] << "," << region._high[1] << "," << region._high[2] << "] with hybrid dim "
        << region._hybridDim  << std::endl;
    }

    _components = _simulation.getEnsemble()->getComponents();
    _domain = domain;
    _particleContainer = particleContainer;

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

        _forceAdapter.init(domain);
    }
}

void AdResS::readXML(XMLfileUnits &xmlconfig) {
    long numRegions = 0;
    XMLfile::Query query = xmlconfig.query("fpregions/region");
    numRegions = query.card();
    if (numRegions == 0) {
        global_log->fatal() << "No AdResS regions specified: Falling back to dynamic region selection. ERROR: not implemented yet!" << std::endl;
        Simulation::exit(668);
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

    for(const auto& region : _fpRegions) {
        global_log->info() << "[AdResS] FPRegion Box from ["
                            << region._low[0] << "," << region._low[1] << "," << region._low[2] << "] to ["
                            << region._high[0] << "," << region._high[1] << "," << region._high[2] << "] with hybrid dim "
                            << region._hybridDim  << std::endl;
    }
}

void AdResS::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                     unsigned long simstep) {

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
        // molecule is in no Hybrid or FP region
        checkMoleculeLOD(*itM, CoarseGrain);
        for(auto& reg : _fpRegions) {
            if(itM->inBox(reg._lowHybrid.data(), reg._highHybrid.data())) {
                checkMoleculeLOD(*itM, Hybrid);
            }
        }
        for(auto& reg : _fpRegions) {
            if(itM->inBox(reg._low.data(), reg._high.data())) {
                checkMoleculeLOD(*itM, FullParticle);
            }
        }
    }
}

void AdResS::siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
    _mesoVals.clear();
    computeForce(true);
    computeForce(false);
    _mesoVals.setInDomain(_domain);
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
    std::array<double, 3> dist = {0,0,0};

    //check all regions
    for(auto& region : _fpRegions) {
        // only molecules within cutoff around the region have interacted with hybrid molecules
        std::array<double, 3> check_low {region._lowHybrid};
        std::array<double, 3> check_high {region._highHybrid};
        for(int d = 0; d < 3; d++) {
            check_low[d] -= cutoff;
            check_high[d] += cutoff;
        }

        auto itOuter = _particleContainer->regionIterator(check_low.data(), check_high.data(), ParticleIterator::ONLY_INNER_AND_BOUNDARY);
        for(; itOuter.isValid(); ++itOuter) {
            Molecule& m1 = *itOuter;

            auto itInner = itOuter;
            ++itInner;
            for(; itInner.isValid(); ++itInner) {
                Molecule& m2 = *itInner;
                mardyn_assert(&m1 != &m2);

                //check if inner is FP or CG -> skip
                if(FPRegion::isInnerPoint(m2.r_arr(), region._low, region._high)) continue;
                if(!FPRegion::isInnerPoint(m2.r_arr(), region._lowHybrid, region._highHybrid)) continue;

                //check distance
                double dd = m1.dist2(m2, dist.data());
                if(dd < cutoff2) {
                    //recompute force and invert it -> last bool param is true
                    _forceAdapter.processPair(m1, m2, dist, MOLECULE_MOLECULE, dd, (dd < LJCutoff2), invert,
                                              _comp_to_res, region);
                }
            }
        }
    }
}

double AdResS::weight(std::array<double, 3> r, FPRegion &region) {
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
