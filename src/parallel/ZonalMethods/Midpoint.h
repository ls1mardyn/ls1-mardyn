/*
 * CommunicationScheme.h
 *
 *  Created on: 14.07.2017
 *      Author: sauermann
 */

#pragma once

#include "Simulation.h"
#include "ZonalMethod.h"

/**
 * This class implements the Midpoint method. Using the Midpoint method the haloRegion is of size r_c/2.
 * Otherwise the halos are identical to the FullShell method.
 * The Midpoint method requires forces to be exchanged.
 */
class Midpoint: public ZonalMethod {
public:
	Midpoint() {
	}
	virtual ~Midpoint() {
	}

	/**
	 * Returns up to 26 halo Regions of the process.
	 * If a process is spanning a whole dimension, then fewer regions can be returned.
	 * The regions indicate, where the processes lie that require halo copies from the current process.
	 * @param initialRegion boundary of the current process
	 * @param cutoffRadius
	 * @return vector of regions
	 */
	virtual std::vector<HaloRegion> getHaloImportForceExportRegions(HaloRegion& initialRegion, double cutoffRadius,
																	bool coversWholeDomain[3],
																	double cellLength[3]) override {
		// the midpoint traversal is cell based, so the halo region has to be cellLength wide.
		return getLeavingExportRegions(initialRegion, cellLength, coversWholeDomain);
	}

	/**
	 * Returns up to 26 halo Regions of the process.
	 * If a process is spanning a whole dimension, then fewer regions can be returned.
	 * The regions indicate, where the processes lie that require halo copies from the current process.
	 * @param initialRegion boundary of the current process
	 * @param cutoffRadius
	 * @return vector of regions
	 */
	virtual std::vector<HaloRegion> getHaloExportForceImportRegions(HaloRegion& initialRegion, double cutoffRadius,
																	bool coversWholeDomain[3],
																	double cellLength[3]) override {
		const std::function<bool(const int[3])> condition = [](const int[3])->bool {
			// no condition for leaving particles.
				return true;
			};
		// the midpoint traversal is cell based, so the halo region has to be cellLength wide.
		return getHaloRegionsConditionalInside(initialRegion, cellLength, coversWholeDomain, condition);
	}
};

