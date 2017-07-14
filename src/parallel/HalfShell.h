/*
 * CommunicationScheme.h
 *
 *  Created on: 14.07.2017
 *      Author: sauermann
 */

#pragma once

#include "CommunicationScheme.h"


class HalfShell: public CommunicationScheme {
public:
	HalfShell(){}
	virtual ~HalfShell(){}

	/**
	 * Returns up to 13 halo Regions of the process. If a process is spanning a whole dimension, then fewer regions can be returned.
	 * @param initialRegion boundary of the current process
	 * @param cutoffRadius
	 * @return vector of regions
	 */
	std::vector<HaloRegion> getHaloRegions(HaloRegion& initialRegion, double cutoffRadius, const bool coversWholeDomain[3]) override;
};

