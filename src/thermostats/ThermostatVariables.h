/*
 * ThermostatVariables.h
 *
 *  Created on: 12 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_THERMOSTATS_THERMOSTATVARIABLES_H_
#define SRC_THERMOSTATS_THERMOSTATVARIABLES_H_

#include <array>

#define USE_THERMOSTAT_VARIABLES 1

class ThermostatVariables {
public:
	ThermostatVariables() { clear(); }
	virtual ~ThermostatVariables() {}

#if 0
	I wanted to disable the copy constructor, but this prevents usage in std containers...
	// disable copy constructor,
	// because forgetting "&" can cause a lot of application programmers a lot of headaches
	ThermostatVariables(const ThermostatVariables& that) = delete;
#endif

	void clear() {
		_numMolecules = 0ul;
		_numRotationalDOF = 0ul;
		_ekinTrans = 0.0;
		_ekinRot = 0.0;
		_betaTrans = 1.0;
		_betaRot = 1.0;
		for (int d = 0; d < 3; ++d) { _totalVelocity[d] = 0.0; }
	}

	void add(const ThermostatVariables& o) {
		_numMolecules += o._numMolecules;
		_numRotationalDOF += o._numRotationalDOF;
		_ekinTrans += o._ekinTrans;
		_ekinRot += o._ekinRot;
	}

	unsigned long _numMolecules, _numRotationalDOF;
	double _ekinTrans, _ekinRot;
	double _betaTrans, _betaRot;
	std::array<double, 3> _totalVelocity;
};

class ThermostatVariablesLocalAndGlobal {
public:
	ThermostatVariablesLocalAndGlobal() { clear(); }
	virtual ~ThermostatVariablesLocalAndGlobal() {}
	void clear() {
		_local.clear();
		_global.clear();
	}
	ThermostatVariables _local, _global;
};

#endif /* SRC_THERMOSTATS_THERMOSTATVARIABLES_H_ */
