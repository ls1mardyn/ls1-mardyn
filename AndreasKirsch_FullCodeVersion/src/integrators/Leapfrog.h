#ifndef LEAPFROG_H_
#define LEAPFROG_H_

#include "integrators/Integrator.h"

//! @brief rotational leapfrog integration scheme
//! @author Martin Buchholz, Martin Bernreuther, et al.
//!
//! For details about the algorithm see David Fincham's paper "Leapfrog rotational algorithms"
//! The integrator responds to the following events:
//!		new timestep -> call preF
//! 	forces calculated -> call postF
//!
//! If you look for the actual leapfrog integration code, go look at preF/postF in Molecule.
//! You ask: Seperation of Concerns? Who cares, fuck it! :)
class Leapfrog : public Integrator {
public:
	//! The constructor
	Leapfrog(double timestepLength);

	//! The destructor
	~Leapfrog();

	//! @brief steps between the force calculation and the end of the time step
	//!
	//! checks whether the current state of the integrator allows that this method is called
	void eventForcesCalculated(ParticleContainer* molCont, Domain* domain);

	//! @brief performs all steps that can be done before new forces are needed
	//!
	//! checks whether the current state of the integrator allows that this method is called
	void eventNewTimestep(ParticleContainer* molCont, Domain* domain);

	virtual void accelerateUniformly(
			ParticleContainer* molCont,
			Domain* domain
	);
	virtual void accelerateInstantaneously(
			ParticleContainer* molCont,
			Domain* domain
	);
	virtual void init1D(
			unsigned zoscillator,
			ParticleContainer* molCont
	);
	virtual void zOscillation(
			unsigned zoscillator,
			ParticleContainer* molCont
	);
};
#endif /*LEAPFROG_H_*/
