/*
 * CommunicationScheme.h
 *
 *  Created on: 14.07.2017
 *      Author: sauermann
 */

#pragma once

#include "FullShell.h"

class Midpoint: public FullShell {
public:
	Midpoint(){}
	virtual ~Midpoint(){}

	/**
	 * Returns up to 26 halo Regions of the process. If a process is spanning a whole dimension, then fewer regions can be returned.
	 * @param initialRegion boundary of the current process
	 * @param cutoffRadius
	 * @return vector of regions
	 */
	std::vector<HaloRegion> getHaloRegions(HaloRegion& initialRegion, double cutoffRadius, const bool coversWholeDomain[3]) override{
		// Same as full shell but with half the cutoff radius
		return FullShell::getHaloRegions(initialRegion, cutoffRadius / 2.0, coversWholeDomain);
	}
};

