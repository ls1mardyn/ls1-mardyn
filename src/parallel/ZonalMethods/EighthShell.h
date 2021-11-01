/*
 * CommunicationScheme.h
 *
 *  Created on: 14.07.2017
 *      Author: sauermann
 */

#pragma once

#include "ZonalMethod.h"

/**
 * This class implements the EighthShell Method. Halo is only imported, if it lies within the eighth shell.
 * The eighth shell is the intrinsic decomposition scheme for any c08 based traversal.
 * The imported region is
 * The EighthShell method requires forces to be exchanged.
 */
class EighthShell: public ZonalMethod {
public:
	EighthShell() = default;
	~EighthShell() override = default;

	std::vector<HaloRegion> getHaloImportForceExportRegions(HaloRegion& initialRegion, double cutoffRadius,
															bool coversWholeDomain[3], double cellLength[3]) override {
		auto condition = [](const int d[3])->bool {
			bool good = true;
			for(unsigned short i = 0; i < 3; ++i){
				good &= d[i] >= 0;
			}
			return good;
		};
		return getHaloRegionsConditional(initialRegion, cutoffRadius, coversWholeDomain, condition);
	}

	std::vector<HaloRegion> getHaloExportForceImportRegions(HaloRegion& initialRegion, double cutoffRadius,
															bool coversWholeDomain[3], double cellLength[3]) override {
		auto condition = [](const int d[3])->bool {
			bool good = true;
			for(unsigned short i = 0; i < 3; ++i){
				good &= d[i] <= 0;
			}
			return good;
		};
		return getHaloRegionsConditionalInside(initialRegion, cutoffRadius, coversWholeDomain, condition);
	}
};

