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

AdResS::AdResS() : _mesoVals(), _fpRegions(), _particleContainer(nullptr), _components(nullptr), _comp_to_res() ,_domain(nullptr) {};

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
    for(auto itM = container->iterator(ParticleIterator::ALL_CELLS); itM.isValid(); ++itM) {
        for(auto& reg : _fpRegions) {
            if(itM->inBox(reg._lowHybrid.data(), reg._highHybrid.data())) {
                // is in FP LOD
                if(itM->inBox(reg._low.data(), reg._high.data())) {
                    checkMoleculeLOD(*itM, FullParticle);
                }
                // is in Hybrid LOD
                else {
                    checkMoleculeLOD(*itM, Hybrid);
                }
                goto end;
            }
        }
        // molecule is in no Hybrid or FP region
        checkMoleculeLOD(*itM, CoarseGrain);

        //molecule has been full handled
        end: continue;
    }
}

void AdResS::siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
    _mesoVals.clear();



    _mesoVals.setInDomain(_domain);
}

void AdResS::checkMoleculeLOD(Molecule &molecule, Resolution targetRes) {
    auto id = molecule.component()->ID();
    //molecule has already correct component
    if(_comp_to_res[id] == targetRes) return;

    int offset = static_cast<int>(targetRes) - static_cast<int>(_comp_to_res[id]);
    molecule.setComponent(&_components->at(id + offset));
}

double AdResS::weight(std::array<double, 3> r, FPRegion &region) {
    std::array<double, 3> axis_vector{ };
    for(int d = 0; d < 3; d++) axis_vector[d] = r[d] - region._center[d];


    return 0;
}
