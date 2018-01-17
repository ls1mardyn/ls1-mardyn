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
	
	void apply(ParticleContainer *moleculeContainer, DomainDecompBase *dode);
	void setGlobalVelocity(double v[3]);
	// handling the scaling factor for the velocity scaling in three dimensions
	void setGlobalBetaTrans(double beta) { _globalBetaTrans = beta; }
	double getGlobalBetaTrans() { return _globalBetaTrans; }
	void setGlobalBetaRot(double beta) { _globalBetaRot = beta; }
	double getGlobalBetaRot() { return _globalBetaRot; }
		
	// handling the scaling factor for the velocity scaling in one dimensions	
	void setGlobalAlphaTrans(double alpha) { _globalAlphaTrans = alpha;}
	double getGlobalAlphaTrans() { return _globalAlphaTrans; }
	
	// handling more than one thermostat
	void enableComponentwise() { _componentwise = true; }
	void disableComponentwise() { _componentwise = false; }
	void setVelocity(int componentId, double v[3]);
	// handling the scaling factor for the velocity scaling in three dimensions
	void setBetaTrans(int componentId, double beta);
	double getBetaTrans(int componentId) { return _componentBetaTrans[componentId]; }
	void setBetaRot(int componentId, double beta);
	double getBetaRot(int componentId) { return _componentBetaRot[componentId]; }
	// handling the scaling factor for the velocity scaling in one dimensions	
	void setAlphaTrans(int componentId, double alpha) { _componentAlphaTrans[componentId] = alpha;}
	double getAlphaTrans(int componentId) { return _componentAlphaTrans[componentId]; }
	
	// calculation of the directed velocities: moving average; size of control volumes: 1*1*z_max OR other
	void calculateDirectedVelocities(ParticleContainer *moleculeContainer, DomainDecompBase *dode);
	// calculation of the directed velocities: simple average; size of control volumes: 1*1*z_max OR other
	void calculateDirectedVelocitiesSimple(ParticleContainer *moleculeContainer, DomainDecompBase *dode);
	// calculation of the directed velocities: simple average; size of control volumes: 1*1*1
	void calculateDirectedVelocitiesSimpleSigma3D(ParticleContainer *moleculeContainer, DomainDecompBase *dode);
	// calculation of the directed velocities: average of neighbour list; size of control volumes: 1*1*z_max OR other
	void calculateDirectedVelocitiesNeighbourList(ParticleContainer *moleculeContainer, DomainDecompBase *dode);

private:
	double _globalBetaTrans;
	double _globalBetaRot;
	double *_globalVelocity;
	double _globalAlphaTrans;
	
	// minimum number of molecules contributing to a valid directed velocity
	unsigned _molDirVelThreshold;

	bool _componentwise;
	std::map<int, double> _componentBetaTrans;
	std::map<int, double> _componentBetaRot;
	std::map<int, double*> _componentVelocity;
	std::map<int, double> _componentAlphaTrans;
	
	// --------------- GENERAL --------------------
	std::map< int, std::map< double, std::map< double, double > > > _directedVelocity[3];
	std::map< int, std::map< double, std::map< double, double > > > _universalDirectedVelocity[3];
	std::map< int, std::map< double, std::map< double, unsigned long > > > _countedMolecules;	
	std::map< int, std::map< double, std::map< double, unsigned long > > > _universalCountedMolecules;
	//corresponding with moving average
	std::map< int, unsigned long > _countedMoleculesAll;
	std::map< int, unsigned long > _universalCountedMoleculesAll;
	std::map< int, double > _directedVelocityAll[3];
	std::map< int, double > _universalDirectedVelocityAll[3];
	std::map< int, std::map< double, std::map< double, unsigned long > > > _averagedN;
	std::map< int, std::map< double, std::map< double, double > > > _averagedVelSum[3];
	std::map< int, std::map< double, std::map< double, double > > > _currentDirectedVel[3];
	std::map< int, std::map< double, std::map< double, std::deque< unsigned long > > > > _universalPriorN;
	std::map< int, std::map< double, std::map< double, std::deque< double > > > > _universalPriorVelSums[3];
	//corresponding with calculateDirectedVelocitiesSimpleSigma3D()
	std::map< int, std::map< double, std::map< double, std::map< double, double > > > > _directedVelocity3D[3];
	std::map< int, std::map< double, std::map< double, std::map< double, double > > > > _universalDirectedVelocity3D[3];
	std::map< int, std::map< double, std::map< double, std::map< double, unsigned long > > > > _countedMolecules3D;
	std::map< int, std::map< double, std::map< double, std::map< double, unsigned long > > > > _universalCountedMolecules3D;
	
	// ------------- CONFINEMENT ------------------
	std::map< int, std::map< unsigned, std::map< unsigned, double > > > _directedVelocityConfinement[3];
	std::map< int, std::map< unsigned, std::map< unsigned, double > > > _universalDirectedVelocityConfinement[3];
	std::map< int, std::map< unsigned, std::map< unsigned, unsigned long > > > _countedMoleculesConfinement;	
	std::map< int, std::map< unsigned, std::map< unsigned, unsigned long > > > _universalCountedMoleculesConfinement;
	//corresponding with moving average
	std::map< int, std::map< unsigned, std::map< unsigned, unsigned long > > > _averagedNConfinement;
	std::map< int, std::map< unsigned, std::map< unsigned, double > > > _averagedVelSumConfinement[3];
	std::map< int, std::map< unsigned, std::map< unsigned, double > > > _currentDirectedVelConfinement[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::deque< unsigned long > > > > _universalPriorNConfinement;
	std::map< int, std::map< unsigned, std::map< unsigned, std::deque< double > > > > _universalPriorVelSumsConfinement[3];
	std::map <unsigned, double> _universalInvProfileUnitConfinement;
	std::map <unsigned, double> _bulkCorner;
	
	// ---------------- SLAB ----------------------
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _directedVelocitySlab[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _universalDirectedVelocitySlab[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, unsigned long > > > > _countedMoleculesSlab;	
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, unsigned long > > > > _universalCountedMoleculesSlab;
	//corresponding with moving average	
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, unsigned long > > > > _averagedNSlab;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _averagedVelSumSlab[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _currentDirectedVelSlab[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::deque< unsigned long > > > > > _universalPriorNSlab;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::deque< double > > > > > _universalPriorVelSumsSlab[3];
	
	// --------------- STRESS ---------------------
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _directedVelocityStress[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _universalDirectedVelocityStress[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, unsigned long > > > > _countedMoleculesStress;	
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, unsigned long > > > > _universalCountedMoleculesStress;
	//corresponding with moving average
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, unsigned long > > > > _averagedNStress;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _averagedVelSumStress[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, double > > > > _currentDirectedVelStress[3];
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::deque< unsigned long > > > > > _universalPriorNStress;
	std::map< int, std::map< unsigned, std::map< unsigned, std::map< unsigned, std::deque< double > > > > > _universalPriorVelSumsStress[3];
	std::map <unsigned, double> _universalInvProfileUnit_Stress;
	
};

#endif /* VELOCITY_SCALING_THERMOSTAT_H */

