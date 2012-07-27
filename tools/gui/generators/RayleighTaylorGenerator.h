/*
 * RayleighTaylorGenerator.h
 *
 *  Created on: June, 2012
 *      Author: nagashim
 */

class Domain;
using namespace std;

#include "MDGenerator.h"
#include "common/ComponentParameters.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompDummy.h"

#ifndef RAYLEIGHTALORGENERATOR_H_
#define RAYLEIGHTALORGENERATOR_H_

class RayleighTaylorGenerator: public MDGenerator {

private://TODO ask Wolfgang why only _temperature is prefixed by _.
	double _temperature;
	int numOfMolecules;
	double L1, L2, L3, epsilon_A, epsilon_B, sigma_A, sigma_B,
		q_A, q_B, m_A, m_B, N, r_cut, delta_t, T, G, h, p_max, skal;
	std::vector<Component> _components;


	double rho;
	double numSphereSizes;

	//double simBoxLength[3]; // we use L1 and L2 in place of simBoxLength[3]

	//! @brief each element is a sphere (vector containing x,y,z and r)
	vector<vector<double> > localClusters;

public:
	RayleighTaylorGenerator();
	virtual ~RayleighTaylorGenerator();

/*	virtual void readPhaseSpaceHeader(Domain* domain, double timestep);

	//! @brief read the phase space components and header information
	unsigned long readPhaseSpace(ParticleContainer* particleContainer,
			std::list<ChemicalPotential>* lmu, Domain* domain,
			DomainDecompBase* domainDecomp);*/

	vector<ParameterCollection*> getParameters();

	//void generatePreview();

	void setParameter(Parameter* p);

	virtual bool validateParameters();

	//! @brief read the phase space components and header information
	//! @param timestep timestep length
	virtual void readPhaseSpaceHeader(Domain* domain, double timestep);

	/**
	 *  @brief read the actual phase space information
	 *  Returns "the highest molecule ID found in the phase space file";
	 *  // todo why? should it be some kind of upper bound for the number of molecules???
	 */
	virtual unsigned long readPhaseSpace(ParticleContainer* particleContainer,
			std::list<ChemicalPotential>* lmu, Domain* domain,
			DomainDecompBase* domainDecomp);

private:

	/**
	 * add a molecule to the container, initializing random velocity, orientation, and so on....
	 */
	void addMolecule(double x, double y, double z, unsigned long id, ParticleContainer* particleContainer);
};

#endif // RayleighTaylorGenerator.h
