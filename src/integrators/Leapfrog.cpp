#include "Leapfrog.h"

#include <map>

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "ensemble/PressureGradient.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/ParticleIterator.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"


using namespace std;
using Log::global_log;

Leapfrog::Leapfrog(double timestepLength) :	Integrator(timestepLength) {
	init();
}

void Leapfrog::init() {
	// set starting state
	_state = STATE_POST_FORCE_CALCULATION;
}

Leapfrog::~Leapfrog() {}

void Leapfrog::readXML(XMLfileUnits& xmlconfig) {
	_timestepLength = 0;
	xmlconfig.getNodeValueReduced("timestep", _timestepLength);
	global_log->info() << "Timestep: " << _timestepLength << endl;
	assert(_timestepLength > 0);
}

void Leapfrog::eventForcesCalculated(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_PRE_FORCE_CALCULATION) {
		transition2to3(molCont, domain);
	}
}

void Leapfrog::eventNewTimestep(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_POST_FORCE_CALCULATION) {
		transition3to1(molCont, domain);
		transition1to2(molCont, domain);
	}
}

void Leapfrog::transition1to2(ParticleContainer* molCont, Domain* /*domain*/) {
	if (this->_state != STATE_NEW_TIMESTEP) {
		global_log->error() << "Leapfrog::transition1to2(...): Wrong state for state transition" << endl;
		return;
	}

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const ParticleIterator begin = molCont->iteratorBegin();
		const ParticleIterator end = molCont->iteratorEnd();

		for (ParticleIterator i = begin; i != end; ++i) {
			i->upd_preF(_timestepLength);
		}
	}

	this->_state = STATE_PRE_FORCE_CALCULATION;
}

void Leapfrog::transition2to3(ParticleContainer* molCont, Domain* domain) {
	if (this->_state != STATE_PRE_FORCE_CALCULATION) {
		global_log->error() << "Leapfrog::transition2to3(...): Wrong state for state transition" << endl;
	}

	/* TODO introduce
	  	class Thermostat {
		unsigned long N, rotDOF;
		double summv2, sumIw2;
	};
	 * but not here, rather in a separate class. The current implementation of the thermostats
	 * (here and in class Domain) is a nightmare.
	 */

	ParticleIterator tM;
	map<int, unsigned long> N;
	map<int, unsigned long> rotDOF;
	map<int, double> summv2;
	map<int, double> sumIw2;
	double dt_half = 0.5 * this->_timestepLength;
	if (domain->severalThermostats()) {
		for (tM = molCont->iteratorBegin(); tM != molCont->iteratorEnd(); ++tM) {
			int cid = tM->componentid();
			int thermostat = domain->getThermostat(cid);
			tM->upd_postF(dt_half, summv2[thermostat], sumIw2[thermostat]);
			N[thermostat]++;
			rotDOF[thermostat] += tM->component()->getRotationalDegreesOfFreedom();
		}
	}
	else {
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			unsigned long Ngt_l = 0;
			unsigned long rotDOFgt_l = 0;
			double summv2gt_l = 0.0;
			double sumIw2gt_l = 0.0;

			const ParticleIterator begin = molCont->iteratorBegin();
			const ParticleIterator end = molCont->iteratorEnd();

			for (ParticleIterator i = begin; i != end; ++i) {
				i->upd_postF(dt_half, summv2gt_l, sumIw2gt_l);
				assert(summv2gt_l >= 0.0);
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
	for (map<int, double>::iterator thermit = summv2.begin(); thermit != summv2.end(); thermit++) {
		domain->setLocalSummv2(thermit->second, thermit->first);
		domain->setLocalSumIw2(sumIw2[thermit->first], thermit->first);
		domain->setLocalNrotDOF(thermit->first, N[thermit->first], rotDOF[thermit->first]);
	}

	this->_state = STATE_POST_FORCE_CALCULATION;

}

void Leapfrog::transition3to1(ParticleContainer* /*molCont*/, Domain* /*domain*/) {
	if (this->_state == STATE_POST_FORCE_CALCULATION) {
		this->_state = STATE_NEW_TIMESTEP;
	}
	else {
		global_log->error() << "Leapfrog::transition3to1(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::accelerateUniformly(ParticleContainer* molCont, Domain* domain) {
	map<unsigned, double>* additionalAcceleration = domain->getPG()->getUAA();
	vector<Component> comp = *(_simulation.getEnsemble()->getComponents());
	vector<Component>::iterator compit;
	map<unsigned, double> componentwiseVelocityDelta[3];
	for (compit = comp.begin(); compit != comp.end(); compit++) {
		unsigned cosetid = domain->getPG()->getComponentSet(compit->ID());
		if (cosetid != 0)
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = _timestepLength * additionalAcceleration[d][cosetid];
		else
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = 0;
	}

	ParticleIterator thismol;
	for (thismol = molCont->iteratorBegin(); thismol != molCont->iteratorEnd(); ++thismol) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(componentwiseVelocityDelta[0][cid],
		              componentwiseVelocityDelta[1][cid],
		              componentwiseVelocityDelta[2][cid]);
	}
}

void Leapfrog::accelerateInstantaneously(ParticleContainer* molCont, Domain* domain) {
	vector<Component> comp = *(_simulation.getEnsemble()->getComponents());
	vector<Component>::iterator compit;
	map<unsigned, double> componentwiseVelocityDelta[3];
	for (compit = comp.begin(); compit != comp.end(); compit++) {
		unsigned cosetid = domain->getPG()->getComponentSet(compit->ID());
		if (cosetid != 0)
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = domain->getPG()->getMissingVelocity(cosetid, d);
		else
			for (unsigned d = 0; d < 3; d++)
				componentwiseVelocityDelta[d][compit->ID()] = 0;
	}

	ParticleIterator thismol;
	for (thismol = molCont->iteratorBegin(); thismol != molCont->iteratorEnd(); ++thismol) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(componentwiseVelocityDelta[0][cid],
		              componentwiseVelocityDelta[1][cid],
		              componentwiseVelocityDelta[2][cid]);
	}
}
