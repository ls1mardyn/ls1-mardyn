//
// Created by alex on 05.04.24.
//

#include "FTH.h"
#include "Simulation.h"
#include "Domain.h"

void FTH::Handler::init(const FTH::Config &config) {
	_config = config;
	if (!config._enableThermodynamicForce) return;
	_config._thermodynamicForceSampleCounter = 0;
}

void FTH::Handler::computeSingleIteration(ParticleContainer &container, const Resolution::FPRegions_t &regions) {
	if(!_config._enableThermodynamicForce) return;
	_config._thermodynamicForceSampleCounter++;
	if(_config._thermodynamicForceSampleCounter % _config._thermodynamicForceSampleGap != 0) return;
	_config._thermodynamicForceSampleCounter = 0;

	if(!_config._createThermodynamicForce) return;
	if(checkConvergence()) {
		Log::global_log->info() << "[AdResS] F_TH has converged." << std::endl;
		writeFinalFTH();
		Simulation::exit(0);
	}
	else updateForce(container, regions);
}

/*****************************************************
 * 					 Grid3DHandler
 * **************************************************/

void FTH::Grid3DHandler::init(const FTH::Config &config) {
	Handler::init(config);
	if (!config._enableThermodynamicForce) return;

	if (_config._fth_file_path.empty()) { // == ""
		// set all FTH values in grid to 0
		for (auto& node : _config._grid->getNodes()) {
			node.data().fth = {0, 0, 0};
		}
	}
	// load existing fth
	else loadGridFTH(*_config._grid, _config._fth_file_path);
}

void FTH::Grid3DHandler::updateForce(ParticleContainer &container, const Resolution::FPRegions_t &regions) {
	if(!_config._enableThermodynamicForce) return;

	// compute gradient of sampled density
	Interpolation::FE::computeGradient<>(*_config._grid,
										 [](node_t& node) { return node.data().density; },
										 [](node_t& node, int dim, double grad) { node.data().grad[dim] = grad; });
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	for (auto & node : _config._grid->getNodes()) {
		auto& node_data = node.data();
		for(int dim = 0; dim < 3; dim++) {
			node_data.fth[dim] -= _config._convergenceFactor * node_data.grad[dim];
		}
	}

	// TODO support multiple regions
	auto& region = regions.at(0);
	auto& low = region._low;
	auto& high = region._high;
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	for (auto & node : _config._grid->getNodes()) {
		const auto& node_pos = node.getPos();
		if (node_pos[0] > low[0] && node_pos[1] > low[1] && node_pos[2] > low[2] &&
			node_pos[0] < high[0] && node_pos[1] < high[1] && node_pos[2] < high[2]) {
			node.data().fth = {0, 0, 0};
		}
	}
}

bool FTH::Grid3DHandler::checkConvergence() {
	if(!_config._enableThermodynamicForce) return false;
	if(_config._convergenceFactor == 0.0) return false;

	auto it = std::max_element(_config._grid->getNodes().begin(), _config._grid->getNodes().end(),
					 [](node_t& n1, node_t& n2) -> bool { return n1.data().density < n2.data().density;});
	double rel_error = std::abs((*it).data().density - _config._rho0) / _config._rho0;
	Log::global_log->info() << "[AdResS] F_TH conv rel error: " << rel_error << std::endl;
	return rel_error <= _config._convergenceThreshold;
}

void FTH::Grid3DHandler::applyForce(ParticleContainer &container, const Resolution::FPRegions_t &regions,
									const Resolution::CompResMap_t &compResMap) {
	if(!_config._enableThermodynamicForce) return;

	std::array<double, 3> low = {0,
								 0,
								 0};
	std::array<double, 3> high= {_simulation.getDomain()->getGlobalLength(0),
								 _simulation.getDomain()->getGlobalLength(1),
								 _simulation.getDomain()->getGlobalLength(2)};
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	for (auto itM = container.regionIterator(std::data(low), std::data(high), ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
		if(compResMap[itM->componentid()] == Resolution::FullParticle) continue;
		if(compResMap[itM->componentid()] == Resolution::CoarseGrain) continue;

		std::array<double, 3> r = itM->r_arr();
		std::array<double, 3> force = interpolateGridFTH(*_config._grid, r);
		itM->Fadd(std::data(force));
	}
}

void FTH::Grid3DHandler::writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain,
								   unsigned long simstep) {
	if(!_config._enableThermodynamicForce) return;
	if(_config._thermodynamicForceSampleCounter != 0) return;

	if (_config._logFTH) writeGridFTH(*_config._grid, "./F_TH_InterpolationFunction", simstep);

	if (!_config._density_sampler->wasSampled()) _config._density_sampler->sampleData(&particleContainer, &domainDecomp, &domain);
	if(_config._logDensities) _config._density_sampler->writeAverage("Density_Grid", simstep);
}

void FTH::Grid3DHandler::writeFinalFTH() {
	if(!_config._enableThermodynamicForce) return;
	writeGridFTH(*_config._grid, "./F_TH_Final");
}

/*****************************************************
 * 					 Grid1DHandler
 * **************************************************/

void FTH::Grid1DHandler::init(const FTH::Config &config) {
	Handler::init(config);
	if (!config._enableThermodynamicForce) return;

	if (_config._fth_file_path.empty()) {
		auto* domain = _simulation.getDomain();
		_thermodynamicForce.n = _bins;
		_thermodynamicForce.begin = 0.0;
		_thermodynamicForce.step_width.resize(_thermodynamicForce.n-1, domain->getGlobalLength(0) / static_cast<double>(_thermodynamicForce.n));
		_thermodynamicForce.gradients.resize(_thermodynamicForce.n, 0.0);
		_thermodynamicForce.function_values.resize(_thermodynamicForce.n, 0.0);
	}
	else {
		_thermodynamicForce.loadTXT(_config._fth_file_path);
	}
}

void FTH::Grid1DHandler::updateForce(ParticleContainer &container, const Resolution::FPRegions_t &regions) {
	if(!_config._enableThermodynamicForce) return;
	auto* density_sampler = _config._density_sampler;

	// compute gradient of sampled density
	auto density = density_sampler->getAverage();
	std::vector<double> d_prime;
	Interpolation::computeGradient(density, d_prime);
	Interpolation::Function d_prime_fun;
	std::vector<double> steps;
	steps.resize(d_prime.size()-1, _thermodynamicForce.step_width[0]);
	Interpolation::computeHermite(0.0, d_prime, steps, d_prime.size(), d_prime_fun);

	for(int i = 0; i < d_prime_fun.n; i++) {
		_thermodynamicForce.function_values[i] -= _config._convergenceFactor * d_prime_fun.function_values[i];
		_thermodynamicForce.gradients[i] -= _config._convergenceFactor * d_prime_fun.gradients[i];
	}

	// TODO support multiple regions
	auto& region = regions.at(0);
	double x_pos = _thermodynamicForce.begin;
	double low = region._low[0];
	double high = region._high[0];
	for(unsigned long i = 0; i < _thermodynamicForce.n; i++) {
		if(x_pos >= low && x_pos <= high){
			_thermodynamicForce.function_values[i] = 0;
			_thermodynamicForce.gradients[i] = 0;
		}
		x_pos += _thermodynamicForce.step_width[i];
	}
}

bool FTH::Grid1DHandler::checkConvergence() {
	if(!_config._enableThermodynamicForce) return false;
	if(_config._convergenceFactor == 0.0) return false;

	auto densities = _config._density_sampler->getAverage();
	auto it = std::max_element(densities.begin(), densities.end());
	double rel_error = std::abs(*it - _config._rho0) / _config._rho0;
	Log::global_log->info() << "[AdResS] F_TH conv rel error: " << rel_error << std::endl;
	return rel_error <= _config._convergenceThreshold;
}

void FTH::Grid1DHandler::applyForce(ParticleContainer &container, const Resolution::FPRegions_t &regions,
									const Resolution::CompResMap_t &compResMap) {
	if(!_config._enableThermodynamicForce) return;

	std::array<double, 3> low = {2*_thermodynamicForce.step_width[0],
								 0,
								 0};
	std::array<double, 3> high= {_simulation.getDomain()->getGlobalLength(0) - 2*_thermodynamicForce.step_width[0],
								 _simulation.getDomain()->getGlobalLength(1),
								 _simulation.getDomain()->getGlobalLength(2)};
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	for (auto itM = container.regionIterator(std::data(low), std::data(high), ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
		if(compResMap[itM->componentid()] == Resolution::FullParticle) continue;
		if(compResMap[itM->componentid()] == Resolution::CoarseGrain) continue;

		double x = itM->r(0);
		double F = computeHermiteAt(x, _thermodynamicForce);
		std::array<double, 3> force = {F, 0.0, 0.0};
		itM->Fadd(std::data(force));
	}
}

void FTH::Grid1DHandler::writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain,
								   unsigned long simstep) {
	if(!_config._enableThermodynamicForce) return;
	if(_config._thermodynamicForceSampleCounter != 0) return;

	std::stringstream stream;
	stream << "./F_TH_InterpolationFunction_" << simstep << ".txt";
	if(_config._logFTH) _thermodynamicForce.writeTXT(stream.str());

	if(!_config._density_sampler->wasSampled()) _config._density_sampler->sampleData(&particleContainer, &domainDecomp, &domain);
	if(_config._logDensities) _config._density_sampler->writeAverage("Density_Grid1D", simstep);
}

void FTH::Grid1DHandler::writeFinalFTH() {
	if(!_config._enableThermodynamicForce) return;
	_thermodynamicForce.writeTXT("./F_TH_Final.txt");
}