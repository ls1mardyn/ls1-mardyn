/*
 * CubicGridGenerator.h
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt
 */

#ifndef EQVGRIDGENERATOR_H_
#define EQVGRIDGENERATOR_H_

#include "common/ComponentParameters.h"
#include "Parameters/ParameterCollection.h"
#include "MDGenerator.h"

#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecompDummy.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "ensemble/PressureGradient.h"
#include "ensemble/GrandCanonical.h"
#include "io/CheckpointWriter.h"
#include "Domain.h"

#include <list>

/**
 * Generate Molecules on a grid according to body centered cubic (bcc) layout,
 * allowing for a density of 68 %.
 */

class CubicGridGenerator: public MDGenerator {

private:
	int _numMolecules;
	double _molarDensity;
	std::vector<Component> _components;
	double _temperature;

	double _simBoxLength[3]; // length of the simulation box

	bool _binaryMixture;
public:
	/**
	 * Constructor
	 */
	CubicGridGenerator();

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
			std::list<ChemicalPotential>* lmu, Domain* domain,
			DomainDecompBase* domainDecomp);

private:
	void calculateSimulationBoxLength();
};

#endif /* EQVGRIDGENERATOR_H_ */
