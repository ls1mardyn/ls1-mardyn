#ifndef INPUTBASE_H_
#define INPUTBASE_H_

#include <string>
#include <list>

class ParticleContainer;

class DomainDecompBase;

class Domain;

class XMLfileUnits;

//! @brief interface for any kind of input class
//!
//! @todo more comment
class InputBase {
public:
	InputBase() = default;

	virtual ~InputBase() = default;

	/** @brief read the phase space components and header information
	 * @param domain  pointer to domain object
	 * @param timestep timestep length
	 */
	virtual void readPhaseSpaceHeader(Domain* domain, double timestep) = 0;

	virtual void readXML(XMLfileUnits& /*xmlconfig*/) {}

	/** @brief read the actual phase space information
	 * @return number of read in molecules added to the particle containers of all processes
	 */
	virtual unsigned long
	readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) = 0;

};

#endif /*INPUTBASE_H_*/
