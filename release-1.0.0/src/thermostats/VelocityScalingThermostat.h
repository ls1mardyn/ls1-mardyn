#ifndef VELOCITY_SCALING_THERMOSTAT_H
#define VELOCITY_SCALING_THERMOSTAT_H

#include "particleContainer/ParticleContainer.h"

#include <map>


class VelocityScalingThermostat {
public:
	VelocityScalingThermostat();
	~VelocityScalingThermostat();
	void setGlobalBetaTrans(double beta) { _globalBetaTrans = beta; }
	double getGlobalBetaTrans() { return _globalBetaTrans; }
	void setGlobalBetaRot(double beta) { _globalBetaRot = beta; }
	double getGlobalBetaRot() { return _globalBetaRot; }
	void setGlobalVelocity(double v[3]);

	void enableComponentwise() { _componentwise = true; }
	void disableComponentwise() { _componentwise = false; }

	void setBetaTrans(int componentId, double beta);
	double getBetaTrans(int componentId) { return _componentBetaTrans[componentId]; }
	void setBetaRot(int componentId, double beta);
	double getBetaRot(int componentId) { return _componentBetaRot[componentId]; }
	void setVelocity(int componentId, double v[3]);
	void apply(ParticleContainer *moleculeContainer);

private:
	double _globalBetaTrans;
	double _globalBetaRot;
	double *_globalVelocity;

	bool _componentwise;
	std::map<int, double> _componentBetaTrans;
	std::map<int, double> _componentBetaRot;
	std::map<int, double*> _componentVelocity;
};

#endif /* VELOCITY_SCALING_THERMOSTAT_H */

