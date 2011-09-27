#ifndef INPUTBASE_H_
#define INPUTBASE_H_

#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase;
class Domain;
class ChemicalPotential;

//! @brief interface for any kind of input class
//!
//! @todo more comment
class InputBase{
public:
	InputBase(){}

	virtual ~InputBase(){}

	//! @brief set the phase space file name
	virtual void setPhaseSpaceFile(std::string filename) = 0;

	//! @brief set the phase space header file name (can be identical to the
	//         phase space file
	virtual void setPhaseSpaceHeaderFile(std::string filename) = 0;

	//! @brief read the phase space components and header information
	//! @param timestep timestep length
	virtual void readPhaseSpaceHeader(Domain* domain, double timestep) = 0;

	/**
	 *  @brief read the actual phase space information
	 *  Returns "the highest molecule ID found in the phase space file";
	 *  // todo why? should it be some kind of upper bound for the number of molecules???
	 */
	virtual unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) = 0;

};

#endif /*INPUTBASE_H_*/
