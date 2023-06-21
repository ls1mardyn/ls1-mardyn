#pragma once

#include <limits>

#include "InputBase.h"
#include "utils/xmlfileUnits.h"

class Component;

/**
 * Class that generates equally distributed particles with exactly the given number of particles per cell.
 * @note: Currently only adds molecules with the first component.
 */
class PerCellGenerator : public InputBase {
public:
	/**
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <generator name="PerCellGenerator">
		<!-- required field: number of molecules per cell -->
		<numMoleculesPerCell>UNSIGNED INT<numMoleculesPerCell/>
		<!-- required field: initial temperature -->
		<initTemperature>DOUBLE<initTemperature/>
		<!-- optional field (default=true): Iff true: Generates two particles per process if numMoleculesPerCell==0.
			 If no particles were generated at all, the simulation crashes.-->
		<generateAtLeastTwoParticles>BOOL</generateAtLeastTwoParticles>
	   </generator>
	   \endcode
	 * @param xmlconfig
	 */
	void readXML(XMLfileUnits &xmlconfig) override;

	void readPhaseSpaceHeader(Domain *domain, double timestep) override{
		// Empty as no header file exists.
	};

	unsigned long readPhaseSpace(ParticleContainer *particleContainer, Domain *domain,
								 DomainDecompBase *domainDecomp) override;

	/**
	 * Fills the container with the specified number of molecules (with given component) per cell.
	 * @param particleContainer
	 * @param component
	 * @param numMoleculesPerCell
	 * @param fillHalo Indicates whether to fill the halo cells as well.
	 */
	static void fillContainer(ParticleContainer *particleContainer, Component *component,
							  unsigned int numMoleculesPerCell, bool fillHalo);

	/**
	 * Generates two particles/molecules at random positions within the boundaries of particleContainer.
	 * @param particleContainer
	 * @param component The component of the particles.
	 */
	static void generateTwoParticles(ParticleContainer *particleContainer, Component *component);

private:
	unsigned int _numMoleculesPerCell{std::numeric_limits<unsigned int>::max()};

	bool _generateAtLeastTwoParticles{true};

	double _initTemperature{0.};
};
