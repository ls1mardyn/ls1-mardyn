/*
 * CubicGridGenerator.h
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt
 */

#ifndef CUBICGRIDGENERATOR_H_
#define CUBICGRIDGENERATOR_H_

#include "common/ComponentParameters.h"
#include "Parameters/ParameterCollection.h"
#include "MDGenerator.h"

#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "io/CheckpointWriter.h"
#include "Domain.h"

#include <list>

/**
 * Generate Molecules on a grid according to body centered cubic (bcc) layout,
 * allowing for a density of 68 %.
 */

class CubicGridGenerator: public MDGenerator {

private:
	// use unsigned long long for BG/P
	unsigned long long int _numMolecules;
	double _molarDensity;
	std::vector<Component>& _components;
	double _temperature;

	double _simBoxLength; // length of the simulation box

	bool _binaryMixture;
public:
	/**
	 * Constructor
	 */
	CubicGridGenerator();

	virtual ~CubicGridGenerator(){};

	/**
	 * Sets the new parameter
	 * @param p the new parameter to be added to the list
	 */
	virtual void setParameter(Parameter* p);

	/**
	 * Creates parameters and puts them into the list of parameters
	 * deletes previously created pa,rameters at first
	 */
	std::vector<ParameterCollection*> getParameters();

	bool validateParameters();

	//! @brief read the phase space components and header information
	//! @param timestep timestep length
	virtual void readPhaseSpaceHeader(Domain* domain, double timestep);

	/**
	 *  @brief read the actual phase space information
	 *  Returns "the highest molecule ID found in the phase space file";
	 *  // todo why? should it be some kind of upper bound for the number of molecules???
	 */
	virtual unsigned long readPhaseSpace(ParticleContainer* particleContainer,
			Domain* domain,
			DomainDecompBase* domainDecomp);

private:
	void calculateSimulationBoxLength();

	/**
	 * add a molecule to the container, initializing random velocity, orientation, and so on....
	 */
	void addMolecule(double x, double y, double z, unsigned long id, ParticleContainer* particleContainer);
};

#endif /* CUBICGRIDGENERATOR_H_ */
