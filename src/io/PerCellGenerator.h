#pragma once

#include "InputBase.h"
#include "utils/xmlfileUnits.h"

class Component;

/**
 * Class that generates equally distributed particles with exactly the given number of particles per cell.
 */
class PerCellGenerator : public InputBase {
public:
	void readPhaseSpaceHeader(Domain *domain, double timestep) override{};

	unsigned long readPhaseSpace(ParticleContainer *particleContainer, Domain *domain,
								 DomainDecompBase *domainDecomp) override;

	/**
	 *
	 * @param particleContainer
	 * @param component
	 * @param numMoleculesPerCell
	 * @param fillHalo
	 */
	static void fillContainer(ParticleContainer *particleContainer, Component *component,
							  unsigned int numMoleculesPerCell, bool fillHalo);

	/**
	 *
	 * @param particleContainer
	 * @param component
	 */
	static void generateTwoParticles(ParticleContainer *particleContainer, Component *component);

	/**
	 *
	 * @param xmlconfig
	 */
	void readXML(XMLfileUnits &xmlconfig) override;

private:

	unsigned int _numMoleculesPerCell{std::numeric_limits<unsigned int>::max()};

	bool _generateAtLeastTwoParticles{true};

	double _initTemperature{0.};
};
