//
// Created by alex on 06.03.24.
//

#include "TemperatureObserver.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"

void TemperatureObserver::readXML(XMLfileUnits &xmlconfig) {
	std::size_t num_regions = 0;
	XMLfile::Query query = xmlconfig.query("ActiveRegions/region");
	num_regions = query.card();
	_region_data.resize(num_regions);

	std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator rIt;
	std::size_t index = 0;
	for (rIt = query.begin(); rIt; rIt++) {
		xmlconfig.changecurrentnode(rIt);
		d3 low{}, high{};
		xmlconfig.getNodeValue("lower/x", low[0]);
		xmlconfig.getNodeValue("lower/y", low[1]);
		xmlconfig.getNodeValue("lower/z", low[2]);

		xmlconfig.getNodeValue("upper/x", high[0]);
		xmlconfig.getNodeValue("upper/y", high[1]);
		xmlconfig.getNodeValue("upper/z", high[2]);

		_region_data[index].first.low = low;
		_region_data[index].first.high = high;
		index++;
	}
	xmlconfig.changecurrentnode(oldpath);

	_max_samples = 100;
	xmlconfig.getNodeValue("samples", _max_samples);
}

void TemperatureObserver::init() {
	for (auto &[region, region_cache]: _region_data) {
		region_cache.history.init(_max_samples);
		region_cache.current_temp = 0.0;
	}
}

void TemperatureObserver::step(ParticleContainer *particleContainer) {
	for (auto &[region, region_cache]: _region_data) {
		region_cache.current_temp = measureTemp(region.low, region.high, particleContainer, region_cache.history);
	}
}

void TemperatureObserver::getRegions(std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> &regions) {
	for (auto &[region, region_cache]: _region_data) {
		regions.emplace_back(region.low, region.high);
	}
}

double TemperatureObserver::getTemperature(std::size_t index) {
	mardyn_assert(index < _region_data.size() && index >= 0);
	return _region_data[index].second.current_temp;
}

double TemperatureObserver::measureTemp(const std::array<double, 3> &low, const std::array<double, 3> &high,
										ParticleContainer *particleContainer, BinData &binData) {
	std::array<double, 3> v_mean = computeMeanVelocityStep(low, high, particleContainer, binData);
	double m = _simulation.getEnsemble()->getComponent(0)->m();
	double v = 0;
	std::size_t n = 0;
	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:v, n)
	#endif
	for (auto it = particleContainer->regionIterator(std::data(low), std::data(high),
													 ParticleIterator::ALL_CELLS); it.isValid(); ++it) {
		for (int d = 0; d < 3; d++) {
			v += std::pow(it->v(d) - v_mean[d], 2);
		}
		n += 1;
	}

	return v * m / (3 * static_cast<double>(n));
}

std::array<double, 3>
TemperatureObserver::computeMeanVelocityStep(const std::array<double, 3> &low, const std::array<double, 3> &high,
											 ParticleContainer *particleContainer, BinData &binData) {
	//compute mean velocity per dim using SMC method, see: Karimian et al. 2011
	std::array<double, 3> v_avg{0};
	double *v_raw = std::data(v_avg);
	std::size_t count = 0;
	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:v_raw[:3], count)
	#endif
	for (auto it = particleContainer->regionIterator(std::data(low), std::data(high),
													 ParticleIterator::ALL_CELLS); it.isValid(); ++it) {
		for (int d = 0; d < 3; d++) {
			v_raw[d] += it->v(d);
		}
		count += 1;
	}
	// setting to 1 in case there were no molecules to not divide by zero; will not change result as v is still zero
	if (count == 0) count = 1;
	binData.mol_counts.insert(count);
	binData.velocities.insert(v_avg);

	//compute CAM average first
	{
		std::size_t steps = binData.mol_counts.size();
		auto d = static_cast<double>(steps <= 1 ? 1 : steps - 1);

		std::array<double, 3> v_cam{0};
		for (std::size_t i = 0; i < steps; i++) {
			const auto &vec = binData.velocities.get(i);
			for (int dim = 0; dim < 3; dim++) {
				v_cam[dim] += vec[dim];
			}
		}
		for (int dim = 0; dim < 3; dim++) {
			v_cam[dim] /= d;
		}

		std::size_t nks = 0;
		for (std::size_t i = 0; i < steps; i++) {
			nks += binData.mol_counts.get(i);
		}
		double nks_d = static_cast<double>(nks) / d;

		for (int dim = 0; dim < 3; dim++) {
			v_cam[dim] /= nks_d;
		}
		binData.v_cam.insert(v_cam);
	}

	//compute SAM over CAM
	{
		std::size_t steps = binData.v_cam.size();
		auto d = static_cast<double>(steps <= 1 ? 1 : steps - 1);

		std::array<double, 3> v_scm{0};
		for (std::size_t i = 0; i < steps; i++) {
			const auto &vec = binData.v_cam.get(i);
			for (int dim = 0; dim < 3; dim++) {
				v_scm[dim] += vec[dim];
			}
		}
		for (int dim = 0; dim < 3; dim++) {
			v_scm[dim] /= d;
		}

		return v_scm;
	}
}
