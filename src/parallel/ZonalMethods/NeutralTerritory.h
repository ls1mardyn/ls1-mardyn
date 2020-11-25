/*
 * NeutralTerritory.h
 *
 *  Created on: Mar 1, 2017
 *      Author: seckler
 */

#pragma once

#include "ZonalMethod.h"

/**
 * This class implements the NeutralTerritory method according to
 * Shaw et. al. (A Fast, Scalable Method ..., 2005)
 *
 */
class NeutralTerritory : public ZonalMethod {
public:
	NeutralTerritory() = default;
	~NeutralTerritory() override = default;

	std::vector<HaloRegion> getHaloImportForceExportRegions(HaloRegion& initialRegion, double cutoffRadius,
															double /*skin*/, bool coversWholeDomain[3],
															double cellLength[3]) override {
		auto condition = [](const int d[3]) -> bool {
			// Determines, whether the region is in the disk.
			// Here, we cannot directly apply the stencil from the NT traversal, as for multiple cells also the
			// direction x=0 and y=-1 needs to be taken into account.
			bool inDisk = (d[2] == 0) && d[0] >= 0;

			// Determines, whether the region is in the tower.
			bool inTower = (d[0] == 0) && (d[1] == 0) && (d[2] != 0);

			// Return true, if region is in the tower or in the disk.
			return inDisk || inTower;
		};
		return getHaloRegionsConditional(initialRegion, cutoffRadius, 0., coversWholeDomain, condition);
	}

	std::vector<HaloRegion> getHaloExportForceImportRegions(HaloRegion& initialRegion, double cutoffRadius,
															double /*skin*/, bool coversWholeDomain[3],
															double cellLength[3]) override {
		auto condition = [](const int d[3]) -> bool {
			// Determines, whether the region is in the disk.
			// Here, we cannot directly apply the stencil from the NT traversal, as for multiple cells also the
			// direction x=0 and y=-1 needs to be taken into account.
			bool inDisk = (d[2] == 0) && d[0] <= 0;
			// Difference to haloimportforceexport: here we have a "<="
			// Determines, whether the region is in the tower.
			bool inTower = (d[0] == 0) && (d[1] == 0) && (d[2] != 0);

			// Return true, if region is in the tower or in the disk.
			return inDisk || inTower;
		};
		return getHaloRegionsConditionalInside(initialRegion, cutoffRadius, 0., coversWholeDomain, condition);
	}
};
