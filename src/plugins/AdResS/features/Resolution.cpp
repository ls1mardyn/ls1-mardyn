//
// Created by alex on 05.04.24.
//

#include "Resolution.h"

#include "Domain.h"
#include "particleContainer/ParticleContainer.h"

void Resolution::Handler::init(Resolution::Config& config) {
	_config = std::move(config);
	for(const auto& region : _config.fpRegions) {
		Log::global_log->debug()
		    << "[AdResS] FPRegion Box from ["
			<< region._low[0] << "," << region._low[1] << "," << region._low[2] << "] to ["
			<< region._high[0] << "," << region._high[1] << "," << region._high[2] << "] with hybrid dims ["
			<< region._hybridDims[0] << "," << region._hybridDims[1] << "," << region._hybridDims[2] << "]" << std::endl;
	}

    // handle components
    if (_config.components->size() != 1) {
        Log::global_log->fatal() << "[AdResS] Only supporting single component simulation at the moment" << std::endl;
        exit(669);
    }

    //=================================
    // create other components
    //=================================
    // Hybrid
    {
        Component comp {};
        comp.setID(1);
        comp.setName("Hybrid");

        // add CG interaction site
        Component& fp_comp = _config.components->at(0);
        std::array<double,3> com{0.0,0.0,0.0};
        const double total_mass = fp_comp.m();
        for (int lj = 0; lj < fp_comp.numLJcenters(); lj++) {
            LJcenter& lj_center = fp_comp.ljcenter(lj);
            for (int dim = 0; dim < 3; dim++) {
                com[dim] += lj_center.m() * lj_center.r()[dim];
            }
        }
        for (int i = 0; i < 3; i++) com[i] /= total_mass;
        LJcenter lj_site {};
        lj_site.setR(0, com[0]); lj_site.setR(1, com[1]); lj_site.setR(2, com[2]);
        lj_site.setM(0);
        std::string site_name = "LJ126";
        lj_site.setName(site_name);
        comp.addLJcenter(lj_site);

        // add all other sites
        for (auto& lj : fp_comp.ljcenters()) {
            LJcenter site_copy = lj;
            comp.addLJcenter(site_copy);
        }
        for (auto& c : fp_comp.charges()) {
            Charge site_copy = c;
            comp.addCharge(site_copy);
        }
        for (auto& d : fp_comp.dipoles()) {
            Dipole site_copy = d;
            comp.addDipole(site_copy);
        }
        for (auto& q : fp_comp.quadrupoles()) {
            Quadrupole site_copy = q;
            comp.addQuadrupole(site_copy);
        }

        comp.setI11(fp_comp.I11());
        comp.setI22(fp_comp.I22());
        comp.setI33(fp_comp.I33());

        _config.components->push_back(comp);
    }
    // CG
    {
        Component comp {};
        comp.setID(2);
        comp.setName("CG");

        Component& fp_comp = _config.components->at(0);
        const double total_mass = fp_comp.m();
        LJcenter lj_site {};
        lj_site.setR(0, 0); lj_site.setR(1, 0); lj_site.setR(2, 0);
        lj_site.setM(total_mass);
        std::string site_name = "LJ126";
        lj_site.setName(site_name);
        comp.addLJcenter(lj_site);

        comp.setI11(fp_comp.I11());
        comp.setI22(fp_comp.I22());
        comp.setI33(fp_comp.I33());

        _config.components->push_back(comp);
    }

    _simulation.getEnsemble()->setComponentLookUpIDs();
    int n = _config.components->size();
    _config.domain->getmixcoeff().resize(n * (n-1), 1.0);
	_config.comp_to_res.resize(_config.components->size());
    if (_config.comp_to_res.size() > 3) {
        Log::global_log->fatal() << "[AdResS] Only supporting single component simulation at the moment" << std::endl;
        exit(669);
    }

	for(Component& comp : *_config.components) {
		unsigned int id = comp.ID();
		if(id % ResolutionCount == FullParticle) _config.comp_to_res[id] = FullParticle;
		if(id % ResolutionCount == Hybrid) _config.comp_to_res[id] = Hybrid;
		if(id % ResolutionCount == CoarseGrain) _config.comp_to_res[id] = CoarseGrain;
	}
}

void Resolution::Handler::checkResolution(ParticleContainer &particleContainer) {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	for(auto itM = particleContainer.iterator(ParticleIterator::ALL_CELLS); itM.isValid(); ++itM) {
		bool stop = false;
		//check if is in any full particle region
		for(auto& reg : _config.fpRegions) {
			if(FPRegion::isInnerPoint(itM->r_arr(), reg._low, reg._high)) {
				checkMoleculeLOD(*itM, FullParticle);
				stop = true;
				break;
			}
		}
		if(stop) continue;

		//check if is in any hybrid region
		for(auto& reg : _config.fpRegions) {
			if(FPRegion::isInnerPoint(itM->r_arr(), reg._lowHybrid, reg._highHybrid)) {
				checkMoleculeLOD(*itM, Hybrid);
				stop = true;
				break;
			}
		}
		if(stop) continue;

		// molecule is in no Hybrid or FP region
		checkMoleculeLOD(*itM, CoarseGrain);
	}
}

void Resolution::Handler::checkMoleculeLOD(Molecule &molecule, ResolutionType targetRes) {
	auto id = molecule.component()->ID();
	//molecule has already correct component
	if(_config.comp_to_res[id] == targetRes) return;

	int offset = static_cast<int>(targetRes) - static_cast<int>(_config.comp_to_res[id]);
	molecule.setComponent(&_config.components->at(id + offset));
}

const Resolution::CompResMap_t &Resolution::Handler::getCompResMap() const {
	return _config.comp_to_res;
}

const Resolution::FPRegions_t &Resolution::Handler::getRegions() const {
	return _config.fpRegions;
}
