#include "thermostats/VelocityScalingThermostat.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"


VelocityScalingThermostat::VelocityScalingThermostat() {
	this->_globalBetaTrans = 1.0;
	this->_globalBetaRot = 1.0;
	this->_componentwise = false;
	this->_useGlobalVelocity = false;
}

VelocityScalingThermostat::~VelocityScalingThermostat() {
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
        this->_useGlobalVelocity = true;
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
	if(_componentwise ) {
		for (auto molecule = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid();
			 ++molecule) {
			int thermostatId;
			double betaTrans = _globalBetaTrans;
			double betaRot = _globalBetaRot;
			int cid = molecule->componentid();
			thermostatId = _simulation.getDomain()->getThermostat(cid);
			betaTrans = _componentBetaTrans[thermostatId];
			betaRot   = _componentBetaRot[thermostatId];

			std::map<int, double*>::iterator vIter;
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
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			double betaTrans = _globalBetaTrans;
			double betaRot = _globalBetaRot;
			global_log->debug() << "Beta rot: " << betaRot << std::endl;
			global_log->debug() << "Beta trans: " << betaTrans << std::endl;

			for (auto i = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); i.isValid(); ++i) {
				if(this->_useGlobalVelocity) {
					i->vsub(_globalVelocity[0], _globalVelocity[1], _globalVelocity[2]);
					i->scale_v(betaTrans);
					i->vadd(_globalVelocity[0], _globalVelocity[1], _globalVelocity[2]);
				}
				else {
					i->scale_v(betaTrans);
				}
				i->scale_D(betaRot);
			}
		} // end pragma omp parallel
	}
}
