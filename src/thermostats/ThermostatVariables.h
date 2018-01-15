/*
 * ThermostatVariables.h
 *
 *  Created on: 12 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_THERMOSTATS_THERMOSTATVARIABLES_H_
#define SRC_THERMOSTATS_THERMOSTATVARIABLES_H_

class LocalThermostatVariables {
public:
	LocalThermostatVariables() : _numMolecules(0ul), _numRotationalDOF(0ul), _ekinTrans(0.0), _ekinRot(0.0) {}
	virtual ~LocalThermostatVariables() {}
	// ThermostatVariables(const LocalThermostatVariables& that) = delete; this prevents "forgetting &, changes not saved" errors, but also usage in std::container-s...

	virtual void clear() { *this = LocalThermostatVariables(); }

	unsigned long _numMolecules, _numRotationalDOF;
	double _ekinTrans, _ekinRot;
};

class GlobalThermostatVariables : public LocalThermostatVariables {
public:
	GlobalThermostatVariables() : LocalThermostatVariables(), _betaTrans(1.0), _betaRot(1.0) {}
	virtual ~GlobalThermostatVariables() {}
	// ThermostatVariables(const GlobalThermostatVariables& that) = delete; this prevents "forgetting &, changes not saved" errors, but also usage in std::container-s...

	virtual void clear() { *this = GlobalThermostatVariables(); }

	double _betaTrans, _betaRot;
};

class LocalAndGlobalThermostatVariables {
public:
	LocalAndGlobalThermostatVariables() : _local(), _global() { }
	virtual ~LocalAndGlobalThermostatVariables() {}

	void clear() { _local.clear(); _global.clear(); }

	LocalThermostatVariables _local;
	GlobalThermostatVariables _global;
};

#endif /* SRC_THERMOSTATS_THERMOSTATVARIABLES_H_ */
