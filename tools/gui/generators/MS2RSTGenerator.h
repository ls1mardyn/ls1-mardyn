/*
 * MS2RSTGenerator.h
 *
 *  Created on: Mar 25, 2011
 *      Author: kovacevt
 */

#ifndef MS2RSTGENERATOR_H_
#define MS2RSTGENERATOR_H_

#include "common/ComponentParameters.h"
#include "Parameters/ParameterCollection.h"
#include "MDGenerator.h"
#include "common/MS2RestartReader.h"

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

class MS2RSTGenerator: public MDGenerator {

private:

	double _molarDensity;
	std::vector<Component> _components;
	double _temperature;

	double _simBoxLength; // length of the simulation box

	// conversion factor from mardyn unit length to angstroem
	double _ms2_to_angstroem;

	std::string _filePath;

	// parameters per component:
	/**
	 * indicate if we have to read in orientation / angular moment / etc.
	 */
	bool _hasRotationalDOF;
	unsigned long long int _numMolecules;

public:
	/**
	 * Constructor
	 */
	MS2RSTGenerator();

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

	void thermostat(ParticleContainer* container);

	/**
	 * add a molecule to the container, initializing random velocity, orientation, and so on....
	 */
	void addMolecule(MS2RestartReader::MoleculeData& ms2mol, ParticleContainer* particleContainer);
};

#endif /* CUBICGRIDGENERATOR_H_ */
