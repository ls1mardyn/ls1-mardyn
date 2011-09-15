/*
 * CubicGridGenerator.h
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt
 */

#ifndef EQVGRIDGENERATOR_H_
#define EQVGRIDGENERATOR_H_

#include "common/DrawableMolecule.h"
#include "QObjects/ScenarioGenerator.h"
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
 * Class CubicGridGenerator implements the generator that generates DrawableMolecules
 * on an equidistant grid
 */

class CubicGridGenerator: public MDGenerator {

private:
	int numMoleculesX, numMoleculesY, numMoleculesZ; // number of mols in each direction
	std::vector<Component> _components;
	double _temperature;
	double dx, dy, dz; // distance between DrawableMolecules
	Position origin; // origin of the first DrawableMolecule

	double simBoxLength[3]; // length of the simulation box

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
};

#endif /* EQVGRIDGENERATOR_H_ */
