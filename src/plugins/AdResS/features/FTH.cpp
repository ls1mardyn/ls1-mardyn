//
// Created by alex on 05.04.24.
//

#include "FTH.h"

FTH::Handler::Handler(const FTH::Config &config) {
	_config._enableThermodynamicForce = config._enableThermodynamicForce;
	if (!config._enableThermodynamicForce) return;
	_config = config;

	auto* domain = _simulation.getDomain();
	_densityProfiler.init(_config._samplingStepSize, domain, _config._smoothingFactor, _config._samplingRadius);
	_config._thermodynamicForceSampleCounter = 0;

	if (!_config._createThermodynamicForce) return;

	_config._thermodynamicForce.n = static_cast<unsigned long>(domain->getGlobalLength(0) / _config._samplingStepSize);
	_config._thermodynamicForce.begin = 0.0;
	_config._thermodynamicForce.step_width.resize(_config._thermodynamicForce.n-1, _config._samplingStepSize);
	_config._thermodynamicForce.gradients.resize(_config._thermodynamicForce.n, 0.0);
	_config._thermodynamicForce.function_values.resize(_config._thermodynamicForce.n, 0.0);
}

void FTH::Handler::computeIteration(ParticleContainer& container, const Resolution::FPRegions_t &regions) {
	if(!_config._enableThermodynamicForce) return;
	//_densityProfiler.sampleDensities(&container, &_simulation.domainDecomposition(), _simulation.getDomain());
	_densityProfiler.computeGMMDensities(&container, &_simulation.domainDecomposition(), _simulation.getDomain());
	//std::vector<double> d{_densityProfiler.getDensitySmoothed(0)};
	Interpolation::Function d_fun = _densityProfiler.getGMMDensity(0);
	Interpolation::Function d_prime_fun;
	Interpolation::computeGradient(d_fun, d_prime_fun);

	//std::vector<double> d_prime;
	//Interpolation::computeGradient(d, d_prime);
	//Interpolation::Function d_prime_fun;
	//std::vector<double> steps;
	//steps.resize(d_prime.size()-1, _config._samplingStepSize);
	//Interpolation::computeHermite(0.0, d_prime, steps, d_prime.size(), d_prime_fun);

	for(int i = 0; i < d_prime_fun.n; i++) {
		_config._thermodynamicForce.function_values[i] -= _config._convergenceFactor * d_prime_fun.function_values[i];
		_config._thermodynamicForce.gradients[i] -= _config._convergenceFactor * d_prime_fun.gradients[i];
	}

	// TODO FIXME!!
	auto& region = regions.at(0);
	double x_pos = _config._thermodynamicForce.begin;
	double low = region._low[0];
	double high = region._high[0];
	for(unsigned long i = 0; i < _config._thermodynamicForce.n; i++) {
		if(x_pos >= low && x_pos <= high){
			_config._thermodynamicForce.function_values[i] = 0;
			_config._thermodynamicForce.gradients[i] = 0;
		}
		x_pos += _config._thermodynamicForce.step_width[i];
	}

	_config._lastGradient = std::move(d_prime_fun);
}

bool FTH::Handler::checkConvergence() {
	if(!_config._enableThermodynamicForce) return false;
	if(_config._lastGradient.function_values.empty()) return false;
	if(_config._convergenceFactor == 0.0) return false;

	std::vector<double> densities{_densityProfiler.getDensitySmoothed(0)};
	auto it = std::max_element(densities.begin(), densities.end());
	double rel_error = std::abs(*it - _config._rho0) / _config._rho0;
	Log::global_log->info() << "[AdResS] F_TH conv rel error: " << rel_error << std::endl;
	return rel_error <= _config._convergenceThreshold;
}

void FTH::Handler::step(ParticleContainer &container, const Resolution::FPRegions_t &regions) {
	if(!_config._enableThermodynamicForce) return;
	_config._thermodynamicForceSampleCounter++;
	_densityProfiler.step(&container);

	if(_config._thermodynamicForceSampleCounter % _config._thermodynamicForceSampleGap != 0) return;
	_config._thermodynamicForceSampleCounter = 0;

	if(!_config._createThermodynamicForce) return;
	if(checkConvergence()) {
		Log::global_log->info() << "[AdResS] F_TH has converged." << std::endl;
		_config._thermodynamicForce.writeXML("./F_TH_InterpolationFunction_Final.xml");
		Simulation::exit(0);
	}
	else computeIteration(container, regions);
}

void FTH::Handler::apply(ParticleContainer &container, const Resolution::FPRegions_t &regions, const Resolution::CompResMap_t& compResMap) {
	if(!_config._enableThermodynamicForce) return;

	std::array<double, 3> low = {2*_config._samplingStepSize,
								 0,
								 0};
	std::array<double, 3> high= {_simulation.getDomain()->getGlobalLength(0) - 2*_config._samplingStepSize,
								 _simulation.getDomain()->getGlobalLength(1),
								 _simulation.getDomain()->getGlobalLength(2)};
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	for (auto itM = container.regionIterator(std::data(low), std::data(high), ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
		if(compResMap[itM->componentid()] == Resolution::FullParticle) continue;
		if(compResMap[itM->componentid()] == Resolution::CoarseGrain) continue;

		double x = itM->r(0);
		double F = computeHermiteAt(x, _config._thermodynamicForce);
		std::array<double, 3> force = {F, 0.0, 0.0};
		itM->Fadd(std::data(force));
	}
}

void FTH::Handler::writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain, unsigned long simstep) {
	if(!_config._enableThermodynamicForce) return;
	if(_config._thermodynamicForceSampleCounter != 0) return;

	std::stringstream stream;
	stream << "./F_TH_InterpolationFunction_" << simstep << ".xml";
	if(_config._logFTH) _config._thermodynamicForce.writeXML(stream.str());

	stream.clear();
	stream = std::stringstream {};
	stream << "./F_TH_Density_" << simstep << ".txt";
	_densityProfiler.sampleDensities(&particleContainer, &domainDecomp, &domain);
	if(_config._logDensities) _densityProfiler.writeDensity(stream.str(), " ", 0, DensityProfile3D::SAMPLE);

	stream.clear();
	stream = std::stringstream {};
	stream << "./F_TH_Density_Smooth_" << simstep << ".txta";
	if(_config._logDensities) _densityProfiler.writeDensity(stream.str(), " ", 0, DensityProfile3D::SMOOTH);

	stream.clear();
	stream = std::stringstream {};
	stream << "./F_TH_Density_GMM_" << simstep << ".xmla";
	_densityProfiler.computeGMMDensities(&particleContainer, &domainDecomp, &domain);
	if(_config._logDensities) _densityProfiler.writeDensity(stream.str(), " ", 0, DensityProfile3D::GMM);

	stream.clear();
	stream = std::stringstream {};
	stream << "./F_TH_Density_FT_" << simstep << ".xmlb";
	_densityProfiler.computeFTDensities(&particleContainer, &domainDecomp, &domain);
	if(_config._logDensities) _densityProfiler.writeDensity(stream.str(), " ", 0, DensityProfile3D::FT);

	stream.clear();
	stream = std::stringstream {};
	stream << "./F_TH_Density_GRID_" << simstep << ".txt";
	if(_config._logDensities) _densityProfiler.writeDensity(stream.str(), " ", 0, DensityProfile3D::GRID);
}

void FTH::Handler::writeFinalFTH() {
	_config._thermodynamicForce.writeXML("./F_TH_Final.xml");
}
