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

	// TODO: This is untested
	std::vector<HaloRegion> getHaloImportForceExportRegions(HaloRegion& initialRegion, double cutoffRadius,
															double /*skin*/, bool coversWholeDomain[3],
															double cellLength[3]) override {
		auto condition = [](const int d[3]) -> bool {
			int pseudoCellIndex = ((d[2] * 2) + d[1]) * 2 + d[2];
			// determines, whether the region is in the disk
			bool inDisk = (d[2] == 0) && pseudoCellIndex > 0;

			// determines, whether the region is in the tower
			bool inTower = (d[0] == 0) && (d[1] == 0) && (d[2] != 0);

			// return true, if region is in the tower or in the disk
			return inDisk || inTower;
		};
		return getHaloRegionsConditional(initialRegion, cutoffRadius, 0., coversWholeDomain, condition);
	}

	// TODO: This is untested
	std::vector<HaloRegion> getHaloExportForceImportRegions(HaloRegion& initialRegion, double cutoffRadius,
															double /*skin*/, bool coversWholeDomain[3],
															double cellLength[3]) override {
		auto condition = [](const int d[3]) -> bool {
			int pseudoCellIndex = ((d[2] * 2) + d[1]) * 2 + d[2];
			// determines, whether the region is in the disk
			bool inDisk = (d[2] == 0) && pseudoCellIndex < 0;
			// difference to haloimportforceexport: here we have a "<"
			// determines, whether the region is in the tower
			bool inTower = (d[0] == 0) && (d[1] == 0) && (d[2] != 0);

			// return true, if region is in the tower or in the disk
			return inDisk || inTower;
		};
		return getHaloRegionsConditional(initialRegion, cutoffRadius, 0., coversWholeDomain, condition);
	}
};
