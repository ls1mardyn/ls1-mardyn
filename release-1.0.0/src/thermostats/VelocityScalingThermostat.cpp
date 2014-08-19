#include "thermostats/VelocityScalingThermostat.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

VelocityScalingThermostat::VelocityScalingThermostat() : _globalBetaTrans(1), _globalBetaRot(1), _globalVelocity(NULL), _componentwise(false) {
	_globalVelocity = new double[3];
}

VelocityScalingThermostat::~VelocityScalingThermostat() {
	delete[] _globalVelocity;
}

void VelocityScalingThermostat::setBetaTrans(int componentId, double beta) {
	_componentBetaTrans[componentId] = beta;
}

void VelocityScalingThermostat::setBetaRot(int componentId, double beta) {
	_componentBetaRot[componentId] = beta;
}

void VelocityScalingThermostat::setGlobalVelocity(double v[3]) {
	for(int d = 0; d < 3; d++) {
		_globalVelocity[d] = v[d];
	}
}
void VelocityScalingThermostat::setVelocity(int componentId, double v[3]) {
	if( _componentVelocity.find(componentId) == _componentVelocity.end() ) {
		_componentVelocity[componentId] = new double[3];
	}
	for(int d = 0; d < 3; d++) {
		_componentVelocity[componentId][d] = v[d];
	}
}

void VelocityScalingThermostat::apply(ParticleContainer *moleculeContainer) {
	Molecule *molecule;
	if(_componentwise ) {
		for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {
			int thermostatId;
			double betaTrans = _globalBetaTrans;
			double betaRot = _globalBetaRot;
			int cid = molecule->componentid();
			thermostatId = _simulation.getDomain()->getThermostat(cid);
			betaTrans = _componentBetaTrans[thermostatId];
			betaRot   = _componentBetaRot[thermostatId];

			map<int, double*>::iterator vIter;
			if( (vIter = _componentVelocity.find(thermostatId)) != _componentVelocity.end()) {
				molecule->scale_v(betaTrans);
			}
			else {
				double *v = vIter->second;
				molecule->vsub(v[0], v[1], v[2]);
				molecule->scale_v(betaTrans);
				molecule->vadd(v[0], v[1], v[2]);
			}
			molecule->scale_D(betaRot);
		}
	}
	else {
		double betaTrans = _globalBetaTrans;
		double betaRot = _globalBetaRot;
		global_log->debug() << "Beta rot: " << betaRot << endl;
		global_log->debug() << "Beta trans: " << betaTrans << endl;
		for (molecule = moleculeContainer->begin(); molecule != moleculeContainer->end(); molecule = moleculeContainer->next()) {
			if(_globalVelocity == NULL) {
				molecule->scale_v(betaTrans);
			}
			else {
				molecule->vsub(_globalVelocity[0], _globalVelocity[1], _globalVelocity[2]);
				molecule->scale_v(betaTrans);
				molecule->vadd(_globalVelocity[0], _globalVelocity[1], _globalVelocity[2]);
			}
			molecule->scale_D(betaRot);
		}
	}
}
