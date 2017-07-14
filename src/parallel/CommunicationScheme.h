/*
 * CommunicationScheme.h
 *
 *  Created on: 14.07.2017
 *      Author: sauermann
 */

#pragma once

#include <vector>

#include "HaloRegion.h"

class CommunicationScheme {
protected:
	CommunicationScheme() {};
public:
	virtual ~CommunicationScheme(){};

	/**
	 * Returns the halo regions of the process. If a process is spanning a whole dimension, fewer regions can be returned.
	 * @param initialRegion boundary of the current process
	 * @param cutoffRadius
	 * @return vector of regions
	 */
	virtual std::vector<HaloRegion> getHaloRegions(HaloRegion& initialRegion, double cutoffRadius, const bool coversWholeDomain[3]) = 0;
};
