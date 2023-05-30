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

double (*AdResS::weight)(std::array<double, 3> r, FPRegion& region) = nullptr;

AdResS::AdResS() : _mesoVals(), _forceAdapter(nullptr), _fpRegions(), _particleContainers(_simulation.getMoleculeContainers()),
                   _components(nullptr), _comp_to_res(), _domain(nullptr) {};

AdResS::~AdResS() = default;

static bool is_init = false;
void AdResS::init(ParticleContainer*, DomainDecompBase *domainDecomp, Domain *domain) {
    if(is_init) return;
    is_init = true;
    global_log->debug() << "[AdResS] Enabled " << std::endl;
    for(const auto& region : _fpRegions) {
        global_log->debug() << "[AdResS] FPRegion Box from ["
        << region._low[0] << "," << region._low[1] << "," << region._low[2] << "] to ["
        << region._high[0] << "," << region._high[1] << "," << region._high[2] << "] with hybrid dims ["
        << region._hybridDims[0] << "," << region._hybridDims[1] << "," << region._hybridDims[2] << "]" << std::endl;
    }

    _components = _simulation.getEnsemble()->getComponents();
    _domain = domain;

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
    static const std::array<std::string, 5> weight_impls {"euclid", "manhattan", "component", "near", "flat"};
    static const std::array<double (*)(std::array<double, 3> r, FPRegion &region), 5> impls
        {weightEuclid, weightManhattan, weightComponent, weightNearest, weightFlat};
    std::string impl = weight_impls[0];
    xmlconfig.getNodeValue("weightImpl", impl);
    unsigned long index = std::distance(weight_impls.begin(), std::find(weight_impls.begin(), weight_impls.end(), impl));
    AdResS::weight = impls[index >= 5 ? 0 : index];
    global_log->info() << "[AdResS] Using weight implementation " << weight_impls[index >= 5 ? 0 : index] << std::endl;

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
    // todo add check that no region overlap even in hybrid considering periodic bounds

    for(const auto& region : _fpRegions) {
        global_log->info() << "[AdResS] FPRegion Box from ["
                            << region._low[0] << "," << region._low[1] << "," << region._low[2] << "] to ["
                            << region._high[0] << "," << region._high[1] << "," << region._high[2] << "] with hybrid dim "
                            << region._hybridDims[0] << "," << region._hybridDims[1] << "," << region._hybridDims[2] << "]" << std::endl;
    }

    //need to set up force handler here
    _forceAdapter = new AdResSForceAdapter(*this);
    _simulation.setParticlePairsHandler(_forceAdapter);
    _domain = _simulation.getDomain();
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
        bool stop = false;

        //check if is in any full particle region
        for(auto& reg : _fpRegions) {
            if(reg.isInnerPointDomain(_domain, FullParticle, itM->r_arr())) {
                checkMoleculeLOD(FullParticle, container, itM);
                stop = true;
                break;
            }
        }
        if(stop) continue;

        //check if is in any hybrid region
        for(auto& reg : _fpRegions) {
            if(reg.isInnerPointDomain(_domain, Hybrid, itM->r_arr())) {
                checkMoleculeLOD(Hybrid, container, itM);
                stop = true;
                break;
            }
        }
        if(stop) continue;

        // molecule is in no Hybrid or FP region
        checkMoleculeLOD(CoarseGrain, container, itM);
    }
    container->updateMoleculeCaches();
}

void AdResS::siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
    if(container != _particleContainers[FullParticle]) return;

    computeForce();
    _mesoVals.setInDomain(_domain);
    _mesoVals.clear();
}

void AdResS::checkMoleculeLOD(Resolution targetRes, ParticleContainer* container, ParticleIterator& it) {
    auto id = it->component()->ID();
    //molecule has already correct component
    if(_comp_to_res[id] == targetRes && container == _particleContainers[targetRes]) return;

    int offset = static_cast<int>(targetRes) - static_cast<int>(_comp_to_res[id]);
    it->setComponent(&_components->at(id + offset));
    Molecule m = *it;
#if defined(_OPENMP)
#pragma omp critical(AdResSCheckMolecule)
    {
        container->deleteMolecule(it, false);
        _particleContainers[targetRes]->addParticle(m, true, true, false);
    }
#endif
}

void AdResS::computeForce() {
    double cutoff = _simulation.getcutoffRadius();
    double cutoff2 = cutoff * cutoff;
    double LJCutoff2 = _simulation.getLJCutoff();
    LJCutoff2 *= LJCutoff2;
    std::array<double,3> globLen{0};
    std::array<double, 3> pc_low{0};
    std::array<double, 3> pc_high{0};
    for(int d = 0; d < 3; d++) {
        globLen[d] = _domain->getGlobalLength(d);
        pc_low[d] = _particleContainers[FullParticle]->getBoundingBoxMin(d);
        pc_high[d] = _particleContainers[FullParticle]->getBoundingBoxMax(d);
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

                    int splitDim = 0;
                    double checkSize = 0;
                    for(int d = 0; d < 3; d++) {
                        if(checkHigh[d] - checkLow[d] > checkSize) {
                            checkSize = checkHigh[d] - checkLow[d];
                            splitDim = d;
                        }
                    }

                    std::array<double,3> off_vec{0,0,0};
                    int cells =std::max(1, static_cast<int>(checkSize/cutoff));
                    off_vec[splitDim] = checkSize / cells;
                    std::vector<RegionParticleIterator> tasks;
                    tasks.resize(3* cells);
                    for(int tID = 0; tID < tasks.size()/3; tID++){
                        std::array<double, 3> start = {}, end = {};
                        for(int d = 0; d < 3; d++) {
                            start[d] = tID * off_vec[d] + checkLow[d];
                            end[d] = (d == splitDim) ? (tID+1) * off_vec[d] + checkLow[d] : checkHigh[d];
                        }
                        tasks[3*tID + FullParticle] = _particleContainers[FullParticle]->regionIterator(start.data(), end.data(), ParticleIterator::ALL_CELLS);
                        tasks[3*tID + Hybrid] = _particleContainers[Hybrid]->regionIterator(start.data(), end.data(), ParticleIterator::ALL_CELLS);
                        tasks[3*tID + CoarseGrain] = _particleContainers[CoarseGrain]->regionIterator(start.data(), end.data(), ParticleIterator::ALL_CELLS);
                    }

                    #if defined(_OPENMP)
                    #pragma omp parallel for
                    #endif
                    for(int tID = 0; tID < tasks.size()/3; tID++){
                        std::array<double, 3> dist = {0,0,0};
                        for(auto itH = tasks[3*tID + Hybrid]; itH.isValid(); ++itH) {
                            Molecule& mH = *itH;
                            for(auto itFP = tasks[3*tID + FullParticle]; itFP.isValid(); ++itFP){
                                Molecule& mFP = *itFP;
                                double dd = mH.dist2(mFP, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mFP, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                            for(auto itCG = tasks[3*tID + CoarseGrain]; itCG.isValid(); ++itCG){
                                Molecule& mCG = *itCG;
                                double dd = mH.dist2(mCG, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mCG, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                        }
                    }

                    #if defined(_OPENMP)
                    #pragma omp parallel for
                    #endif
                    for(int tID = 0; tID < (tasks.size()-3)/3; tID+=2){
                        std::array<double, 3> dist = {0,0,0};
                        for(auto itH = tasks[3*tID + Hybrid]; itH.isValid(); ++itH) {
                            Molecule& mH = *itH;
                            for(auto itFP = tasks[3*tID + FullParticle +3]; itFP.isValid(); ++itFP){
                                Molecule& mFP = *itFP;
                                double dd = mH.dist2(mFP, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mFP, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                            for(auto itCG = tasks[3*tID + CoarseGrain +3]; itCG.isValid(); ++itCG){
                                Molecule& mCG = *itCG;
                                double dd = mH.dist2(mCG, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mCG, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                        }
                        for(auto itH = tasks[3*tID + Hybrid +3]; itH.isValid(); ++itH) {
                            Molecule& mH = *itH;
                            for(auto itFP = tasks[3*tID + FullParticle]; itFP.isValid(); ++itFP){
                                Molecule& mFP = *itFP;
                                double dd = mH.dist2(mFP, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mFP, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                            for(auto itCG = tasks[3*tID + CoarseGrain]; itCG.isValid(); ++itCG){
                                Molecule& mCG = *itCG;
                                double dd = mH.dist2(mCG, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mCG, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                        }
                    }

                    #if defined(_OPENMP)
                    #pragma omp parallel for
                    #endif
                    for(int tID = 1; tID <= (tasks.size()-6)/3; tID+=2){
                        std::array<double, 3> dist = {0,0,0};
                        for(auto itH = tasks[3*tID + Hybrid]; itH.isValid(); ++itH) {
                            Molecule& mH = *itH;
                            for(auto itFP = tasks[3*tID + FullParticle +3]; itFP.isValid(); ++itFP){
                                Molecule& mFP = *itFP;
                                double dd = mH.dist2(mFP, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mFP, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                            for(auto itCG = tasks[3*tID + CoarseGrain +3]; itCG.isValid(); ++itCG){
                                Molecule& mCG = *itCG;
                                double dd = mH.dist2(mCG, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mCG, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                        }
                        for(auto itH = tasks[3*tID + Hybrid +3]; itH.isValid(); ++itH) {
                            Molecule& mH = *itH;
                            for(auto itFP = tasks[3*tID + FullParticle]; itFP.isValid(); ++itFP){
                                Molecule& mFP = *itFP;
                                double dd = mH.dist2(mFP, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mFP, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                            for(auto itCG = tasks[3*tID + CoarseGrain]; itCG.isValid(); ++itCG){
                                Molecule& mCG = *itCG;
                                double dd = mH.dist2(mCG, dist.data());
                                if(dd < cutoff2) {
                                    _forceAdapter->processPair(mH, mCG, dist.data(), MOLECULE_MOLECULE, dd, (dd < LJCutoff2), _comp_to_res, false, region);
                                }
                            }
                        }
                    }

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
