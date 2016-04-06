/*
 * AqueousNaClGenerator.h
 *
 *  Created on: Oct 14, 2012
 *      Author: eckhardw
 */

#ifndef AQUEOUSNACLGENERATOR_H_
#define AQUEOUSNACLGENERATOR_H_

#include "common/ComponentParameters.h"
#include "Parameters/ParameterCollection.h"
#include "MDGenerator.h"

#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "ensemble/PressureGradient.h"
#include "ensemble/GrandCanonical.h"
#include "io/CheckpointWriter.h"
#include "Domain.h"

#include <list>

/**
 * Generate a 1 molar aqueous NaCl solution of water. This is a scenario to check
 * calculation of long-range interactions and has been used in
 *
 * Tironi, Sperb, Smith and van Gunsteren:
 * A generalized reaction field method for molecular dynamics simulations
 * J. Chem. Phys. 102, 5451 (1995), doi: 10.1063/1.469273
 *
 * The generation of the initial configuration is done differently:
 * - the number of ions is calculated from the total number of molecules
 * - the ions are randomly placed on the lattice
 * - the water molecules are placed on the lattice if there's no ion already.
 */
class AqueousNaClGenerator: public MDGenerator {

private:

	unsigned long long int _numMolecules;

	std::vector<Component> _components;
	double _temperature;

	double _simBoxLength; // length of the simulation box
	double _molarDensity;
public:
	/**
	 * Constructor
	 */
	AqueousNaClGenerator();

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

	struct Ion {
		int cid;
		int position[3];
	};

	void calculateSimulationBoxLength();

	bool isNearIon(Ion ions[], int ni) const;

	int getCID(Ion ions[], int numIons, int x, int y, int z) const;

	/**
	 * add a molecule to the container, initializing random velocity, orientation, and so on....
	 */
	void addMolecule(double x, double y, double z, unsigned long id, int cid, ParticleContainer* particleContainer);
};

#endif /* AQUEOUSNACLGENERATOR_H_ */
