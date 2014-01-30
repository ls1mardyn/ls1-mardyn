#include "Leapfrog.h"

#include <map>

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "ensemble/PressureGradient.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
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

void Leapfrog::transition1to2(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_NEW_TIMESTEP) {
		Molecule* tempMolecule;
		double vcorr = 2. - 1. / domain->getGlobalBetaTrans();
		double Dcorr = 2. - 1. / domain->getGlobalBetaRot();
		for (tempMolecule = molCont->begin(); tempMolecule != molCont->end(); tempMolecule = molCont->next()) {
			tempMolecule->upd_preF(_timestepLength, vcorr, Dcorr);
		}

		this->_state = STATE_PRE_FORCE_CALCULATION;
	}
	else {
		global_log->error() << "Leapfrog::transition1to2(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::transition2to3(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_PRE_FORCE_CALCULATION) {
		Molecule* tM;
		map<int, unsigned long> N;
		map<int, unsigned long> rotDOF;
		map<int, double> summv2;
		map<int, double> sumIw2;
		double dt_half = 0.5 * this->_timestepLength;
		if (domain->severalThermostats()) {
			for (tM = molCont->begin(); tM != molCont->end(); tM = molCont->next()) {
				int cid = tM->componentid();
				int thermostat = domain->getThermostat(cid);
				tM->upd_postF(dt_half, summv2[thermostat], sumIw2[thermostat]);
				N[thermostat]++;
				rotDOF[thermostat] += tM->component()->getRotationalDegreesOfFreedom();
			}
		}
		else {
			unsigned long Ngt = 0;
			unsigned long rotDOFgt = 0;
			double summv2gt = 0.0;
			double sumIw2gt = 0.0;
			for (tM = molCont->begin(); tM != molCont->end(); tM = molCont->next()) {
				tM->upd_postF(dt_half, summv2gt, sumIw2gt);
				assert(summv2gt >= 0.0);
				Ngt++;
				rotDOFgt += tM->component()->getRotationalDegreesOfFreedom();
			}
			N[0] = Ngt;
			rotDOF[0] = rotDOFgt;
			summv2[0] = summv2gt;
			sumIw2[0] = sumIw2gt;
		}
		for (map<int, double>::iterator thermit = summv2.begin(); thermit != summv2.end(); thermit++) {
			domain->setLocalSummv2(thermit->second, thermit->first);
			domain->setLocalSumIw2(sumIw2[thermit->first], thermit->first);
			domain->setLocalNrotDOF(thermit->first, N[thermit->first], rotDOF[thermit->first]);
		}

		this->_state = STATE_POST_FORCE_CALCULATION;
	}
	else {
		global_log->error() << "Leapfrog::transition2to3(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::transition3to1(ParticleContainer* molCont, Domain* domain) {
	if (this->_state == STATE_POST_FORCE_CALCULATION) {
		this->_state = STATE_NEW_TIMESTEP;
	}
	else {
		global_log->error() << "Leapfrog::transition3to1(...): Wrong state for state transition" << endl;
	}
}

void Leapfrog::accelerateUniformly(ParticleContainer* molCont, Domain* domain) {
	map<unsigned, double>* additionalAcceleration = domain->getPG()->getUAA();
	vector<Component> comp = *(_simulation.getEnsemble()->components());
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

	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(componentwiseVelocityDelta[0][cid],
		              componentwiseVelocityDelta[1][cid],
		              componentwiseVelocityDelta[2][cid]);
	}
}

void Leapfrog::accelerateInstantaneously(ParticleContainer* molCont, Domain* domain) {
	vector<Component> comp = *(_simulation.getEnsemble()->components());
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

	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next()) {
		unsigned cid = thismol->componentid();
		assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
		thismol->vadd(componentwiseVelocityDelta[0][cid],
		              componentwiseVelocityDelta[1][cid],
		              componentwiseVelocityDelta[2][cid]);
	}
}

void Leapfrog::init1D(unsigned zoscillator, ParticleContainer* molCont) {
	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
		if (!(thismol->id() % zoscillator) && thismol->numTersoff()) thismol->setXY();
}

void Leapfrog::zOscillation(unsigned zoscillator, ParticleContainer* molCont) {
	Molecule* thismol;
	for (thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
		/* TODO: use cid instead of complicated id + tersoff */
		if (!(thismol->id() % zoscillator) && thismol->numTersoff()) thismol->resetXY();
}
