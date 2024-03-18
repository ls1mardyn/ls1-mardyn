/*
 * RayleighTaylorGenerator.h
 *
 * For the problem description, please see
 * Numerical Simulation in Molecular Dynamics (NSMD) page 281 - page 284
 * Initial parameter values are take from Table 7.1 on page 284
 *
 * For reference values, see NSMD page 96
 *
 *  Created on: June, 2012
 *      Author: nagashim
 */

class Domain;

#include "MDGenerator.h"
#include "common/ComponentParameters.h"
#include "parallel/DomainDecompBase.h"

#ifndef RAYLEIGHTALORGENERATOR_H_
#define RAYLEIGHTALORGENERATOR_H_

class RayleighTaylorGenerator: public MDGenerator {

private:
//	double _temperature;
	int  _n_1, _n_2, _n_3;
	double _L1, _L2, _L3, _epsilon_A, _epsilon_B, _sigma_A, _sigma_B,
		_q_A, _q_B, _m_A, _m_B, _T;
	std::vector<Component>& _components;

	double numSphereSizes;

public:
/*
	static const double angstroem_2_atomicUnitLength;

	static const double unitMass_2_mardyn; // Mardyn calculates with 1/1000 u as base unit.

	static const double debye_2_mardyn;

	static const double buckingham_2_mardyn;

	static const double unitCharge_2_mardyn;*/

	// tildes are defined on page 282 of NSMD
	static constexpr double atomic_mass_unit_u 	= 1.6605655e-27;//[kg]
	static constexpr double m_tilde				= 1. * 1.6605655e-27; //[kg]

	static constexpr double elementary_charge	= 1.6021892e-19;//[C]
	static constexpr double q_tilde				= 1. * 1.6021892e-19;//[C]

	static constexpr double sigma_tilde 		= 2.22; //[ngstr m]
	static constexpr double epsilon_tilde 		= 1.04710e-21; // [J]

	/**
	 * Constructor
	 */
	RayleighTaylorGenerator();
	virtual ~RayleighTaylorGenerator();

	std::vector<ParameterCollection*> getParameters();

	void setParameter(Parameter* p);

	virtual bool validateParameters();

	//! @brief read the phase space components and header information
	//! @param timestep timestep length
	virtual void readPhaseSpaceHeader(Domain* domain, double timestep);


	//
	//  @brief read the actual phase space information
	//  Returns "the highest molecule ID found in the phase space file";
	//  // todo why? should it be some kind of upper bound for the number of molecules???
	//
	virtual unsigned long readPhaseSpace(ParticleContainer* particleContainer,
			Domain* domain,
			DomainDecompBase* domainDecomp);

private:
	//
	// add a molecule to the container, initializing random velocity, orientation, and so on....
	//
	void addMolecule(double x, double y, double z, unsigned long id,
			int componentType, ParticleContainer* particleContainer);
};

#endif // RayleighTaylorGenerator.h

