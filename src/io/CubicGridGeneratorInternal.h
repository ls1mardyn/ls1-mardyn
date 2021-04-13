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

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {
	}
	unsigned long readPhaseSpace(ParticleContainer *particleContainer, Domain *domain, DomainDecompBase *domainDecomp);

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

	/**
	 * determine Nx, Ny, Nz, s.t.
	 * Nx * Ny * Nz * 2 = targetTotalNumMols
	 * and
	 * Nx : Ny : Nz = boxLength[0] : boxLength[1] : boxLength[2]
	 */
	std::array<unsigned long, 3> determineMolsPerDimension(unsigned long targetTotalNumMols, std::array<double, 3> boxLength) const;

};
