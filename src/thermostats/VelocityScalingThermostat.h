#ifndef VELOCITY_SCALING_THERMOSTAT_H
#define VELOCITY_SCALING_THERMOSTAT_H

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"

#include <map>
#include <deque>

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
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > > _directedVelocitySlab;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > > _universalDirectedVelocitySlab;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > > _directedVelocityStress;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > > _universalDirectedVelocityStress;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _directedVelocityConfinement;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _universalDirectedVelocityConfinement;
	//std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _directedVelocityConfinementTest;
	//std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > >_universalDirectedVelocityConfinementTest;
	std::map< int, std::map< unsigned, std::map< unsigned, double > > > _universalDirectedVelocityTest;
	std::map< int, std::map< double, std::map< double, unsigned long > > > _countedMolecules;	
	std::map< int, std::map< double, std::map< double, unsigned long > > > _universalCountedMolecules;
	
	std::map< int, unsigned long > _countedMoleculesAll;
	std::map< int, unsigned long > _universalCountedMoleculesAll;
	std::map< int, double > _directedVelocityAll[3];
	std::map< int, double > _universalDirectedVelocityAll[3];
	
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< double, unsigned long > > > > _countedMoleculesSlab;	
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< double, unsigned long > > > > _universalCountedMoleculesSlab;	
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< double, unsigned long > > > > _countedMoleculesStress;	
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< double, unsigned long > > > > _universalCountedMoleculesStress;	
	std::map< int, std::map< unsigned, std::map< unsigned, unsigned long > > > _countedMoleculesConfinement;	
	std::map< int, std::map< unsigned, std::map< unsigned, unsigned long > > > _universalCountedMoleculesConfinement;
	//std::map< int, std::map< unsigned, std::map< unsigned, unsigned long > > > _countedMoleculesConfinementTest;
	//std::map< int, std::map< unsigned, std::map< unsigned, unsigned long > > > _universalCountedMoleculesConfinementTest;
	std::map< int, std::map< double, std::map< double, unsigned long > > > _averagedN;
	std::map< int, std::map< double, std::map< double, std::map< unsigned, double > > > > _averagedVelSum;
	std::map< int, std::map< double, std::map< double, std::map< unsigned, double > > > > _currentDirectedVel;
	std::map< int, std::map< double, std::map< double, std::deque< unsigned long > > > > _universalPriorN;
	std::map< int, std::map< double, std::map< double, std::deque< double > > > > _universalPriorVelSums[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, unsigned long > > > > _averagedNSlab;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > > _averagedVelSumSlab;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > > _currentDirectedVelSlab;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::deque< unsigned long > > > > > _universalPriorNSlab;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::deque< double > > > > > _universalPriorVelSumsSlab[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, unsigned long > > > > _averagedNStress;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > > _averagedVelSumStress;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > > _currentDirectedVelStress;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::deque< unsigned long > > > > > _universalPriorNStress;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::deque< double > > > > > _universalPriorVelSumsStress[3];
	std::map< int, std::map< unsigned, std::map< unsigned, unsigned long > > > _averagedNConfinement;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _averagedVelSumConfinement;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _currentDirectedVelConfinement;
	std::map< int, std::map< unsigned, std::map< unsigned, std::deque< unsigned long > > > > _universalPriorNConfinement;
	std::map< int, std::map< unsigned, std::map< unsigned, std::deque< double > > > > _universalPriorVelSumsConfinement[3];
		
	std::map <unsigned, double> _universalInvProfileUnit_Stress;
	std::map <unsigned, double> _bulkCorner;
	std::map <unsigned, double> _universalInvProfileUnitConfinement;
};

#endif /* VELOCITY_SCALING_THERMOSTAT_H */

