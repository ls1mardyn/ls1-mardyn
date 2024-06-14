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

	_config.comp_to_res.resize(_config.components->size());
	for(Component& comp : *_config.components) {
		unsigned int id = comp.ID();
		if(id % ResolutionCount == FullParticle) {
			if(comp.getName().rfind("FP_", 0) == 0) _config.comp_to_res[id] = FullParticle;
			else {
				Log::global_log->fatal() << "[AdResS] Component with id=" << id+1 << " does not have FP_ prefix in name."
									<< "Uncertain if this component should be Full Particle representation." << std::endl;
				exit(669);
			}
		}
		if(id % ResolutionCount == Hybrid) {
			if(comp.getName().rfind("H_", 0) == 0) _config.comp_to_res[id] = Hybrid;
			else {
				Log::global_log->fatal() << "[AdResS] Component with id=" << id+1 << " does not have H_ prefix in name."
									<< "Uncertain if this component should be Hybrid representation." << std::endl;
				exit(669);
			}
		}
		if(id % ResolutionCount == CoarseGrain) {
			if(comp.getName().rfind("CG_", 0) == 0) _config.comp_to_res[id] = CoarseGrain;
			else {
				Log::global_log->fatal() << "[AdResS] Component with id=" << id+1 << " does not have CG_ prefix in name."
									<< "Uncertain if this component should be Coarse Grain representation." << std::endl;
				exit(669);
			}
		}
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
