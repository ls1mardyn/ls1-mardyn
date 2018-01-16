/*
 * FullShell.h
 *
 *  Created on: Oct 10, 2016
 *      Author: seckler
 */

#pragma once

#include "ZonalMethod.h"


class FullShell: public ZonalMethod {
public:
	FullShell(){}
	virtual ~FullShell(){}

	/**
	 * Returns up to 26 halo Regions of the process.
	 * If a process is spanning a whole dimension, then fewer regions can be returned.
	 * The regions indicate, where the processes lie that require halo copies from the current process.
	 * @param initialRegion boundary of the current process
	 * @param cutoffRadius
	 * @return vector of regions
	 */
	virtual std::vector<HaloRegion> getHaloImportForceExportRegions(HaloRegion& initialRegion, double cutoffRadius,
			bool coversWholeDomain[3]) override {
		return getLeavingExportRegions(initialRegion, cutoffRadius, coversWholeDomain);
	}

	virtual std::vector<HaloRegion> getHaloExportForceImportRegions(HaloRegion& initialRegion, double cutoffRadius,
				bool coversWholeDomain[3]) override {
		const std::function<bool(const int[3])> condition = [](const int[3])->bool {
			// no condition for leaving particles.
				return true;
			};
		return getHaloRegionsConditionalInside(initialRegion, cutoffRadius, coversWholeDomain, condition);
	}
};

