#ifndef INPUTBASE_H_
#define INPUTBASE_H_

#include <string>
#include <list>

class ParticleContainer;
class DomainDecompBase;
class Domain;
class ChemicalPotential;
class XMLfileUnits;

//! @brief interface for any kind of input class
//!
//! @todo more comment
class InputBase{
public:
	InputBase(){}

	virtual ~InputBase(){}

	//! @brief read the phase space components and header information
	//! @param timestep timestep length
	virtual void readPhaseSpaceHeader(Domain* domain, double timestep) = 0;

	virtual void readXML(XMLfileUnits& /*xmlconfig*/) {}

	/**
	 *  @brief read the actual phase space information
	 *  Returns "the highest molecule ID found in the phase space file";
	 *  // todo why? should it be some kind of upper bound for the number of molecules???
	 *  // WE: good question. For the scenario generator however, useful for the visualization...
	 */
	virtual unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) = 0;

};

#endif /*INPUTBASE_H_*/
