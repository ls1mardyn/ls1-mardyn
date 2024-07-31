/*
 * MDGenerator.h
 *
 * @Date: 08.06.2011
 * @Author: eckhardw
 */

#ifndef MDGENERATOR_H_
#define MDGENERATOR_H_

#include "common/MardynConfiguration.h"

#include "Generators/Generator.h"
#include "io/InputBase.h"
#include "molecules/Component.h"
#include "utils/Logger.h"

#include <iostream>
#include <list>
#include <vector>
#include <memory>

class Domain;
class DomainDecompBase;
class ParticleContainer;
class ParameterCollection;
class Parameter;

/**
 * This is the base class for all molecule generators and offers functionality
 * which is common to all of them, like writing Mardyn input files or creating
 * drawable molecules.
 *
 * @note if macro MARDYN is defined, this class is built for use within mardyn,
 *       otherwise it is built for use from the ScenarioGenerator. Both modes are
 *       exclusive.
 */
class MDGenerator: public Generator, public InputBase {

public:

	// multiplicative factors for unit conversions
	// the factors are based on the file atomic_units.txt
	static const double angstroem_2_atomicUnitLength;

	static const double unitMass_2_mardyn; // Mardyn calculates with 1/1000 u as base unit.

	static const double debye_2_mardyn;

	static const double buckingham_2_mardyn;

	static const double unitCharge_2_mardyn;

	static const double abogadro_constant;

	static const double boltzmann_constant_kB;

	/**
	 * femto-seconds to mardyn
	 */
	static const double fs_2_mardyn;

	/**
	 * convert particle density mol/liter to particles / \f$a_0^3\f$
	 *
	 * 1 mol / l = 1 mol / 10^{-3} m^3
	 *           = 6.022141 * 10^{23} / 10^{-3} * (1.889726 * 10^{10} a_0)^3
	 *           = 6.022141 * 10^{23} / 6.748333 * 10^{27} a_0^3
	 *           = 8.923894 * 10^{-5} / a_0^3
	 */
	static const double molPerL_2_mardyn;

	static const double kelvin_2_mardyn;

protected:
	std::shared_ptr<Log::Logger> _logger;

	MDGenerator(std::string name);

	virtual ~MDGenerator() {};

	/**
	 * determine the velocity according to the temperature.
	 */
	std::vector<double> getRandomVelocity(double temperature) const;

	/**
	 * Determine if the given position is inside the domain.
	 */
	bool isInsideDomain(Domain* domain, double position[3]);

public:

	virtual void setLogger(std::shared_ptr<Log::Logger> logger);

	//! NOP
	void setPhaseSpaceFile(std::string /*filename*/) {}

	//! NOP
	void setPhaseSpaceHeaderFile(std::string /*filename*/) {}

	/**
	 * Generates DrawableMolecules and saves them in the list
	 */
	virtual void generatePreview();

	/**
	 * Writes output file(s)
	 */
	virtual void generateOutput(const std::string& directory);

	void createSampleObject() const;

	/***********************************************/
	/*** Methods to be implemented by subclasses ***/
	/***********************************************/

	/**
	 * Creates parameters and puts them into the list of parameters.
	 * The caller is responsible for destructing the objects returned!
	 *
	 * @return generator parameters
	 */
	virtual std::vector<ParameterCollection*> getParameters() = 0;

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
	virtual unsigned long readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) = 0;

	/**
	 * Remove the system of momentum from particle container.
	 */
	static void removeMomentum(ParticleContainer* particleContainer, const std::vector<Component>& components);

protected:

	MardynConfiguration _configuration;

	/**
	 * create a random number between a and b (inclusive)
	 */
	double randdouble(double a, double b) const {
		return a + rand() * (b - a) / (RAND_MAX);
	}

	void getOrientation(int base, int delta, double orientation[4]);

};

#endif /* MDGENERATOR_H_ */
