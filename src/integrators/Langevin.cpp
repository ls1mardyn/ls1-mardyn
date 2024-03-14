//
// Created by alex on 05.03.24.
//

#include "Langevin.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"
#include "utils/xmlfileUnits.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "Domain.h"

#include <random>

void Langevin::readXML(XMLfileUnits &xmlconfig) {
	_checkFailed = false;
	_timestepLength = 0;
	xmlconfig.getNodeValueReduced("timestep", _timestepLength);
	Log::global_log->info() << "Timestep: " << _timestepLength << std::endl;
	mardyn_assert(_timestepLength > 0);

	_gamma = 0;
	xmlconfig.getNodeValue("friction", _gamma);
	mardyn_assert(_gamma > 0);

	std::size_t num_regions = 0;
	XMLfile::Query query = xmlconfig.query("ActiveRegions/region");
	num_regions = query.card();
	_stochastic_regions.resize(num_regions);

	std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator rIt;
	std::size_t index = 0;
	for(rIt = query.begin(); rIt; rIt++) {
		xmlconfig.changecurrentnode(rIt);
		d3 low {}, high{};
		xmlconfig.getNodeValue("lower/x", low[0]);
		xmlconfig.getNodeValue("lower/y", low[1]);
		xmlconfig.getNodeValue("lower/z", low[2]);

		xmlconfig.getNodeValue("upper/x", high[0]);
		xmlconfig.getNodeValue("upper/y", high[1]);
		xmlconfig.getNodeValue("upper/z", high[2]);

		_stochastic_regions[index].low = low;
		_stochastic_regions[index].high = high;
		index++;
	}
	xmlconfig.changecurrentnode(oldpath);
}

void Langevin::init() {
	_dt_half = _timestepLength / 2;

	if (_simulation.getTemperatureObserver() == nullptr && !_checkFailed) {
		_checkFailed = true;
		return;
	} // we will check again later

	if (_simulation.getTemperatureObserver() == nullptr && _checkFailed) {
		Log::global_log->warning() << "Langevin Integrator used, but no Temperature Observer defined as thermostat."
								   << " Behaviour is undefined."
								   << " Use Temperature Observer to disable other thermostat functionality!" << std::endl;
	}
}

void Langevin::eventForcesCalculated(ParticleContainer *particleContainer, Domain *domain) {
	addLangevinContribution(particleContainer);

	std::map<int, unsigned long> N;
	std::map<int, unsigned long> rotDOF;
	std::map<int, double> summv2;
	std::map<int, double> sumIw2;

	if (domain->severalThermostats()) {
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			std::map<int, unsigned long> N_l;
			std::map<int, unsigned long> rotDOF_l;
			std::map<int, double> summv2_l;
			std::map<int, double> sumIw2_l;

			for (auto tM = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
				int cid = tM->componentid();
				int thermostat = domain->getThermostat(cid);
				tM->upd_postF(_dt_half, summv2_l[thermostat], sumIw2_l[thermostat]);
				N_l[thermostat]++;
				rotDOF_l[thermostat] += tM->component()->getRotationalDegreesOfFreedom();
			}

			#if defined(_OPENMP)
			#pragma omp critical (thermostat)
			#endif
			{
				for (auto &it: N_l) N[it.first] += it.second;
				for (auto &it: rotDOF_l) rotDOF[it.first] += it.second;
				for (auto &it: summv2_l) summv2[it.first] += it.second;
				for (auto &it: sumIw2_l) sumIw2[it.first] += it.second;
			}
		}
	} else {
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			unsigned long Ngt_l = 0;
			unsigned long rotDOFgt_l = 0;
			double summv2gt_l = 0.0;
			double sumIw2gt_l = 0.0;

			for (auto i = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); i.isValid(); ++i) {
				i->upd_postF(_dt_half, summv2gt_l, sumIw2gt_l);
				mardyn_assert(summv2gt_l >= 0.0);
				Ngt_l++;
				rotDOFgt_l += i->component()->getRotationalDegreesOfFreedom();
			}

			#if defined(_OPENMP)
			#pragma omp critical (thermostat)
			#endif
			{
				N[0] += Ngt_l;
				rotDOF[0] += rotDOFgt_l;
				summv2[0] += summv2gt_l;
				sumIw2[0] += sumIw2gt_l;
			}
		} // end pragma omp parallel
	}
	for (auto &thermit: summv2) {
		domain->setLocalSummv2(thermit.second, thermit.first);
		domain->setLocalSumIw2(sumIw2[thermit.first], thermit.first);
		domain->setLocalNrotDOF(thermit.first, N[thermit.first], rotDOF[thermit.first]);
	}
}

void Langevin::eventNewTimestep(ParticleContainer *particleContainer, Domain *domain) {
	addLangevinContribution(particleContainer);

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		for (auto i = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); i.isValid(); ++i) {
			i->upd_preF(_timestepLength);
		}
	}


}

Langevin::d3 Langevin::sampleRandomForce(double m, double T) {
	static thread_local std::mt19937 generator; // keeping default seed for reproducibility
	std::normal_distribution normal{0.0, 1.0};
	d3 r_vec{};

	// additional factor 2 here
	// according to book about Langevin integration not needed, but according to Langevin's equations of motion it does exist
	// by using it we actually reach the target temperature
	double scale = std::sqrt(2 * _timestepLength * T * _gamma / m);
	for (int d = 0; d < 3; d++) {
		r_vec[d] = scale * normal(generator);
	}

	return r_vec;
}

void Langevin::addLangevinContribution(ParticleContainer *particleContainer) {
	for (std::size_t index = 0; index < _stochastic_regions.size(); index++) {
		auto &region = _stochastic_regions[index];
		double T = _simulation.getEnsemble()->T();
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		for (auto it = particleContainer->regionIterator(std::data(region.low),
														 std::data(region.high),
														 ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
			d3 v_old = it->v_arr();
			d3 v_rand = sampleRandomForce(it->mass(), T);
			d3 v_delta{};

			for (int d = 0; d < 3; d++) {
				v_delta[d] = -_dt_half * _gamma * v_old[d] + v_rand[d];
			}

			it->vadd(v_delta[0], v_delta[1], v_delta[2]);
		}
	}
}
