/*
 * FullShell.h
 *
 *  Created on: Oct 10, 2016
 *      Author: seckler
 */

#pragma once

#include <vector>

#include "HaloRegion.h"


class FullShell {
public:
	FullShell();
	virtual ~FullShell();
	/**
	 * Returns up to 26 halo Regions of the process. If a process is spanning a whole dimension, then fewer regions can be returned.
	 * @param initialRegion boundary of the current process
	 * @param cutoffRadius
	 * @return vector of regions
	 */
	std::vector<HaloRegion> getHaloRegions(HaloRegion initialRegion, double cutoffRadius, const bool coversWholeDomain[3]);
};

