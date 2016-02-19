#ifndef VELOCITY_SCALING_THERMOSTAT_H
#define VELOCITY_SCALING_THERMOSTAT_H

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"

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
	
	void setGlobalAlphaTrans(double alpha) { _globalAlphaTrans = alpha;}
	double getGlobalAlphaTrans() { return _globalAlphaTrans; }
	void setAlphaTrans(int componentId, double alpha) { _componentAlphaTrans[componentId] = alpha;}
	double getAlphaTrans(int componentId) { return _componentAlphaTrans[componentId]; }

	void enableComponentwise() { _componentwise = true; }
	void disableComponentwise() { _componentwise = false; }

	void setBetaTrans(int componentId, double beta);
	double getBetaTrans(int componentId) { return _componentBetaTrans[componentId]; }
	void setBetaRot(int componentId, double beta);
	double getBetaRot(int componentId) { return _componentBetaRot[componentId]; }
	void setVelocity(int componentId, double v[3]);
	void apply(ParticleContainer *moleculeContainer, DomainDecompBase *dode);
	void calculateDirectedVelocities(ParticleContainer *moleculeContainer, DomainDecompBase *dode);


private:
	double _globalBetaTrans;
	double _globalBetaRot;
	double *_globalVelocity;
	double _globalAlphaTrans;

	bool _componentwise;
	std::map<int, double> _componentBetaTrans;
	std::map<int, double> _componentBetaRot;
	std::map<int, double*> _componentVelocity;
	std::map<int, double> _componentAlphaTrans;
	
	std::map< int, std::map< double, std::map< double, std::map< unsigned, double > > > > _directedVelocity;
	std::map< int, std::map< double, std::map< double, std::map< unsigned, double > > > > _universalDirectedVelocity;
	std::map< int, std::map< unsigned, std::map< unsigned, double > > > _universalDirectedVelocityTest;
	std::map< int, std::map< double, std::map< double, unsigned long> > > _countedMolecules;	
	std::map< int, std::map< double, std::map< double, unsigned long > > > _universalCountedMolecules;	
	std::map <unsigned, double> _universalInvProfileUnit_Stress;
	std::map <unsigned, double> _bulkCorner;
	std::map <unsigned, double> _universalInvProfileUnitConfinement;
};

#endif /* VELOCITY_SCALING_THERMOSTAT_H */

