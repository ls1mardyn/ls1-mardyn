#ifndef LEAPFROG_H_
#define LEAPFROG_H_

#include "integrators/Integrator.h"

//! @brief rotational leapfrog integration scheme
//! @author Martin Buchholz, Martin Bernreuther, et al.
//!
//! For details about the algorithm see David Fincham's paper "Leapfrog rotational algorithms"
//! This Leapfrog integrator is implemented as a deterministic finite automaton (DFA):
//! - A state of the DFA corresponds to a coarse step in the integration loop.
//!   The following states are used:
//!   - state 3: "starting state". The simulation has either just started, or the last steps
//!     of a time step have just been done. In state 3, the leapfrog integrator waits for the start
//!     of the next time time step
//!   - state 1: a new time step has just begun, the integrator is ready to start integrating
//!     (transition to state 2)
//!   - state 2: the first part is finished, now forces are needed to continue the integration
//! - consequently, the transitions do the following things:
//!   - 3 to 1: Nothing to do
//!   - 1 to 2: calculate "preF"
//!   - 2 to 3: calculate "postF"
//! - For the transition to another state, a condition has to be fulfilled
//!     There are two possibilities: The integrator has to wait for external information
//!     e.g. the new forces. In this case, the integrator is informed, that new forces
//!     have been calculated. The automaton can then do all necessary computations
//!     that have to be done do get from the current state to the next state
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

private:

	//! state in which the integrator is
	int _state;

	//! @brief calculate new positions and the first velocity halfstep
	//!
	//! This method also checks whether the state is 1. If so, the calculations are done and
	//! the state is set to 2, otherwise, an error is printed
	void transition1to2(ParticleContainer* molCont, Domain* domain);

	//! @brief calculates the second halfstep for the velocity and angular momentum (postF)
	//!
	//! This method also checks whether the state is 2. If so, the calculations are done and
	//! the state is set to 3, otherwise, an error is printed
	void transition2to3(ParticleContainer* molCont, Domain* domain);
	//! @brief checks whether the state is 3. If so, the state is set to 3, otherwise, an error is printed.
	void transition3to1(ParticleContainer* molCont, Domain* domain);

};
#endif /*LEAPFROG_H_*/
