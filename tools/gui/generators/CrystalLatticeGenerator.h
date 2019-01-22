/*
 * CrystalLatticeGenerator.h
 *
 *  Created on: Oct 13, 2012
 *      Author: Wolfgang Eckhardt
 */

#ifndef CRYSTALLATTICEGENERATOR_H_
#define CRYSTALLATTICEGENERATOR_H_

#include "common/ComponentParameters.h"
#include "Parameters/ParameterCollection.h"
#include "MDGenerator.h"

#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "ensemble/GrandCanonical.h"
#include "io/CheckpointWriter.h"
#include "Domain.h"

#include <list>

/**
 * Generate Molecules in a crystal lattice, reprecenting a structure as it is
 * found for NaCl. In this scenario, the electrostatic potential energy can be calculated
 * analytically via the Madelung constant.
 *
 * This is a static scenario, i.e. all forces should euqal 0.
 */

class CrystalLatticeGenerator: public MDGenerator {

private:

	unsigned long long int _numMoleculesPerDim;
	// lattice spacing
	double _h;
	double _charge;
	std::vector<Component>& _components;

	double _simBoxLength; // length of the simulation box

public:
	/**
	 * Constructor
	 */
	CrystalLatticeGenerator();

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

	/**
	 * add a molecule to the container, initializing random velocity, orientation, and so on....
	 */
	void addMolecule(double x, double y, double z, unsigned long id, unsigned cid, ParticleContainer* particleContainer);
};

#endif /* CRYSTALLATTICEGENERATOR_H_ */
