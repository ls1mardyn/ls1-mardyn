#ifndef OUTPUTBASE_H_
#define OUTPUTBASE_H_

#include "ensemble/GrandCanonical.h"
#include <list>

class ParticleContainer;
class DomainDecompBase; 
class Domain;

/**
 * TODO: cleanup all classes implementing this interface
 *       make this an abstract class which deals with all things which are
 *       common to all OutputComponents like writeFrequncy, baseFilename, and so on...
 */

//! @brief interface for any kind of output class
//! @author Martin Bernreuther, Martin Buchholz, et al. (2009)
//!
//! There are lots of different things that possibly might want
//! to be written out, e.g. thermodynimic values, graphical information, 
//! time measurements, ... \n
//! In many cases, the output happens regularly every time step. To keep
//! the other classes as clear as possible, the output routines (especially
//! extensive routines writing special graphics formats for visualization)
//! shouldn't be placed there. That's why this interface was introduced.
//! For a given output task, this interface has to be implemented. This
//! implementation will be called OutputPlugin in the follwing. 
//! 
//! There are three main methodes
//! - initOutput: will be called once in the beginning
//! - doOutput: will be called each time step
//! - finishOutput: will be called at the end
//! To each of this methods, a pointer to the particle Container, to the domain 
//! decomposition and to the domain will be passed. This should cover most of
//! the data that might be necessary for the output.
//! 
//! Of course, several OutputPlugins will be needed in some cases. So the
//! idea is, that the class which controls the simulation (see Simulation.h)
//! has a list of OutputPlugins. At the beginning of the simulation, 
//! for each element in that list the method initOutput is called. 
//! The same will happen in each time step with the method doOutput and at 
//! the end of the simulation with the method finishOutput.
class OutputBase {
public:
	//! @brief Subclasses should use their constructur to pass parameters (e.g. filenames)
	OutputBase(){}

	virtual ~OutputBase(){}

	//! @brief will be called at the beginning of the simulation
	//!
	//! Some OutputPlugins will need some initial things to be done before
	//! the output can start, e.g. opening some files. This method will
	//! be called once at the beginning of the simulation (see Simulation.cpp)
	virtual void initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) = 0;

	//! @brief will be called in each time step
	//!
	//! Most of the times, the output should either be done every time step or at least
	//! every n-th time step. Therefore, this method has an additional parameter simstep,
	//! allowing to do a output depending on the current simulation time step. This method
	//! will be called once every time step during the simulation (see Simulation.cpp)
	virtual void doOutput(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep,
			std::list<ChemicalPotential>* lmu
	) = 0;

	//! @brief will be called at the end of the simulation
	//!
	//! Some OutputPlugins will need to do some things at the end of the simulation,
	//! e.g. closing some files. This method will
	//! be called once at the end of the simulation (see Simulation.cpp)
	virtual void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) = 0;
};

#endif /*OUTPUTBASE_H_*/
