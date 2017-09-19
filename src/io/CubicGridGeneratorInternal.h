/*
 * CubicGridGeneratorInternal.h
 *
 *  Created on: Jun 9, 2017
 *      Author: seckler
 */

#pragma once
#include <vector>
#include "InputBase.h"
#include "molecules/Component.h"
#include "utils/Random.h"

class ParticleContainer;
class ChemicalPotential;
class Domain;
class DomainDecompBase;


class CubicGridGeneratorInternal: public InputBase {
private:
	// use unsigned long long for BG/P
	unsigned long long int _numMolecules;
	bool _binaryMixture;

public:
	CubicGridGeneratorInternal();
	virtual ~CubicGridGeneratorInternal() {
	}

	void setPhaseSpaceFile(std::string /*filename*/) {
	}
	void setPhaseSpaceHeaderFile(std::string /*filename*/) {
	}

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {
	}
	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu,
			Domain* domain, DomainDecompBase* domainDecomp);

	/** @brief Read in XML configuration for MkTcTSGenerator and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	 <generator name="CubicGridGenerator">
	 TODO: more
	 </generator>
	 \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);
private:
	bool addMolecule(double x, double y, double z, unsigned long id, ParticleContainer* particleContainer);
	void removeMomentum(ParticleContainer* particleContainer, const std::vector<Component>& components);
	/**
	 * create a random number between a and b (inclusive)
	 */
	double randdouble(double a, double b) {
		return _RNG.uniformRandInRange(a, b);
	}
	void getOrientation(int base, int delta, double orientation[4]);

	/**
	 * determine the velocity according to the temperature.
	 */
	std::vector<double> getRandomVelocity(double temperature);

	Random _RNG;
};
