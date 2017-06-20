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
#include "molecules/Molecule.h"

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
	void bufferMolecule(double x, double y, double z, unsigned long id, ParticleContainer* particleContainer);
	void insertMoleculesInContainer(ParticleContainer* particleContainer);
	void removeMomentum(ParticleContainer* particleContainer, const std::vector<Component>& components);
	/**
	 * create a random number between a and b (inclusive)
	 */
	double randdouble(double a, double b) const {
		return a + rand() * (b - a) / (RAND_MAX);
	}
	void getOrientation(int base, int delta, double orientation[4]);

	/**
	 * determine the velocity according to the temperature.
	 */
	std::vector<double> getRandomVelocity(double temperature) const;

	static const int _MOLECULE_BUFFER_SIZE = 1000;
	std::vector<Molecule> _moleculeBuffer;
};
