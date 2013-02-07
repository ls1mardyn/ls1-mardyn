#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

class ParticleContainer;
class Domain;

//! @brief Update velocities and positions.
//! @author Martin Buchholz, Martin Bernreuther, et al.
//!
//! The basic idea of the integrator is, that it calculates from
//! given values at one timestep the corresponding values in the next timestep.
//! A PRECONDITION for each integrator is, that all values that the integrator
//! reads before it writes them have to be available at the same point in
//! time (usually t=0). E.g. positions and velocities are always available
//! at the same starting time as they are read from the input file, but
//! the forces don't stand in the input file, so they have to be calculated
//! for t=0 before the integrator starts working.
//!
//! The idea of this interface is that different integrators can be used,
//! but the user of the main application shouldn't have to care which
//! integrator is used. This of course means, that there can't be a
//! method called "second_velocity_half_step_of_leapfrog_integrator".
//! For different integrators, different steps at different steps of
//! the program are needed. So it isn't feasible to have methods which
//! performs a certain step of an integrator. The "User" even might not know,
//! which information the integrator needs to calculate a certain value.
//! So the idea is, the the user (surrounding method) informs the integrator
//! whenever some new information (e.g. new Forces on the molecules) is
//! available. The integrator should then know what can be done with
//! this new information (e.g. calculate new accelaration)
class Integrator {
public:
	//! The constructor
	Integrator();

	//! The destructor
	virtual ~Integrator();

	//! @brief informs the integrator about available forces, who then continues integration
	//!
	//! An Integrator can't calculate the force on molecules. But the forces are needed to solve
	//! the equations of motion. This method informs the integrator that the forces on all molecules
	//! have been calculated (by some other module), so the integrator should continue it's work
	//! up to a point where all values are available at the time step end.
	//! @param moleculeContainer containes the molecules for which the equations of motion shall be solved
	//! @param domain needed because some macroscopic values (Thermostat) might influence the integrator
	virtual void eventForcesCalculated(ParticleContainer* moleculeContainer, Domain* domain) = 0;

	//! @brief inform the integrator that the integration should continue
	//!
	//! Basically, the only point where an integrator has to interrupt it's work is when
	//! it needs forces to continue, which is not necessarily between two time steps. But
	//! the simulation program usually wants to do something (output) between two time steps
	//! (this has to be done as only then all values (x, v,...) are avaiable for the same point in time).
	//! So the integrator interrupts it's work at the end of one time step and resumes work
	//! with the call of this method.
	//! @param moleculeContainer containes the molecules for which the equations of motion shall be solved
	//! @param domain needed because some macroscopic values (Thermostat) might influence the integrator
	virtual void eventNewTimestep(ParticleContainer* moleculeContainer, Domain* domain) = 0;

	//! set the time between two time steps
	void setTimestepLength(double dt) {
		_timestepLength = dt;
	}

	//! get the time between two time steps
	double getTimestepLength() {
		return _timestepLength;
	}

	virtual void accelerateUniformly(
			ParticleContainer* molCont,
			Domain* domain
	) = 0;
	virtual void accelerateInstantaneously(
			ParticleContainer* molCont,
			Domain* domain
	) = 0;
	virtual void init1D(
			unsigned zoscillator,
			ParticleContainer* molCont
	) = 0;
	virtual void zOscillation(
			unsigned zoscillator,
			ParticleContainer* molCont
	) = 0;

protected:

	//! time between time step n and time step (n+1)
	double _timestepLength;

};
#endif /*INTEGRATOR_H_*/
