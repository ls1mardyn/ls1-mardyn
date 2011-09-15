/*
 * MDGenerator.h
 *
 * @Date: 08.06.2011
 * @Author: eckhardw
 */

#ifndef MDGENERATOR_H_
#define MDGENERATOR_H_

#include "Generators/Generator.h"
#include "common/MardynConfiguration.h"
#include "io/InputBase.h"
#include "molecules/Component.h"
#include "utils/Logger.h"

#include <iostream>
#include <list>

class Domain;
class DomainDecompBase;
class ParticleContainer;
class ChemicalPotential;


/**
 * This is the base class for all molecule generators and offers functionality
 * which is common to all of them, like writing Mardyn input files or creating
 * drawable molecules.
 *
 * @note if macro MARDYN is defined, this class is built for use within mardyn,
 *       otherwise it is built for use from the ScenarioGenerator. Both modes are
 *       exclusive.
 */
class MDGenerator: public Generator, InputBase {

protected:
	Log::Logger* _logger;
	bool _deleteLogger;

	MDGenerator(std::string name);

	virtual ~MDGenerator();

	/**
	 * determine the velocity according to the temperature.
	 */
	std::vector<double> getRandomVelocity(double temperature) const;

	/**
	 * Determine if the given position is inside the domain.
	 */
	bool isInsideDomain(Domain* domain, double position[3]);

public:

	void setLogger(Log::Logger* logger);

	//! NOP
	void setPhaseSpaceFile(std::string filename) {}

	//! NOP
	void setPhaseSpaceHeaderFile(std::string filename) {}

	/**
	 * Generates DrawableMolecules and saves them in the list
	 */
	virtual void generatePreview();

	/**
	 * Writes output file(s)
	 */
	virtual void generateOutput(const std::string& directory);

	const Object* getSampleObject() const;

	/***********************************************/
	/*** Methods to be implemented by subclasses ***/
	/***********************************************/

	/**
	 * Creates parameters and puts them into the list of parameters.
	 * The caller is responsible for destructing the objects returned!
	 *
	 * @return generator parameters
	 */
	virtual vector<ParameterCollection*> getParameters() = 0;

	/**
	* Sets a new parameter and adds it to the list
	*/
	virtual	void setParameter(Parameter* p) = 0;

	/**
	 * Validates if parameters are ok
	 */
	virtual bool validateParameters() = 0;

	//! @brief read the phase space components and header information
	//! @param timestep timestep length
	virtual void readPhaseSpaceHeader(Domain* domain, double timestep) = 0;

	/**
	 *  @brief read the actual phase space information
	 *  Returns "the highest molecule ID found in the phase space file";
	 *  // todo why? should it be some kind of upper bound for the number of molecules???
	 */
	virtual unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) = 0;

protected:

	MardynConfiguration _configuration;

private:

	double randdouble(double a, double b) const {
		return a + rand() * (b - a) / (RAND_MAX);
	}
};

#endif /* MDGENERATOR_H_ */
