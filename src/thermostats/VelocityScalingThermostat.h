#ifndef VELOCITY_SCALING_THERMOSTAT_H
#define VELOCITY_SCALING_THERMOSTAT_H

#include "particleContainer/ParticleContainer.h"

#include <map>


class VelocityScalingThermostat {
public:
	VelocityScalingThermostat();
	~VelocityScalingThermostat();
	void setGlobalBetaTrans(double beta) { _globalBetaTrans = beta; }
	void setGlobalBetaRot(double beta) { _globalBetaRot = beta; }
	void setGlobalVelocity(double v[3]);

	void enableComponentwise() { _componentwise = true; }
	void disableComponentwise() { _componentwise = false; }

	void setBetaTrans(int componentId, double beta);
	void setBetaRot(int componentId, double beta);
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