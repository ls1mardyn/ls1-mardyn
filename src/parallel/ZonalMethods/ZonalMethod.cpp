/*
 * ZonalMethod.cpp
 *
 *  Created on: Feb 28, 2017
 *      Author: seckler
 */

#include "ZonalMethod.h"
#include "Simulation.h"

ZonalMethod::ZonalMethod() {

}

ZonalMethod::~ZonalMethod() {

}

std::vector<HaloRegion> ZonalMethod::getLeavingExportRegions(HaloRegion& initialRegion, double cutoffRadius,
		bool coversWholeDomain[3]) {
	const std::function<bool(const int[3])> condition = [](const int[3])->bool {
		// no condition for leaving particles.
		return true;
	};
	return getHaloRegionsConditional(initialRegion, cutoffRadius, coversWholeDomain, condition);
}



// protected, used for child classes
std::vector<HaloRegion> ZonalMethod::getHaloRegionsConditional(HaloRegion& initialRegion, double cutoffRadius,
			bool coversWholeDomain[3], const std::function<bool(const int[3])>& condition){
	std::vector<HaloRegion> regions;
		int d[3];
		for (d[0] = -1; d[0] <= 1; d[0]++) {
			for (d[1] = -1; d[1] <= 1; d[1]++) {
				for (d[2] = -1; d[2] <= 1; d[2]++) {
					if ((d[0] || d[1] || d[2]) == 0)  // if all are 0 (false), then continue
						continue;

					// we don't include anything, that does not conform with the condition
					if (!condition(d)) {
						continue;
					}

					HaloRegion tmp = initialRegion;
					for (unsigned int dimension = 0; dimension < 3; dimension++) {
						if (d[dimension] == 0) {
							tmp.rmin[dimension] = initialRegion.rmin[dimension];
							tmp.rmax[dimension] = initialRegion.rmax[dimension];
						} else if (d[dimension] == -1) { // LOWER
							tmp.rmin[dimension] = initialRegion.rmin[dimension] - cutoffRadius;
							tmp.rmax[dimension] = initialRegion.rmin[dimension];
						} else { //d[dimension==1 - UPPER
							tmp.rmin[dimension] = initialRegion.rmax[dimension];
							tmp.rmax[dimension] = initialRegion.rmax[dimension] + cutoffRadius;
						}
					}
					tmp.offset[0] = d[0];
					tmp.offset[1] = d[1];
					tmp.offset[2] = d[2];

					tmp.width = cutoffRadius;
					regions.push_back(tmp);
				}
			}
		}
		return regions;
}

// protected, used for child classes
std::vector<HaloRegion> ZonalMethod::getHaloRegionsConditionalInside(HaloRegion& initialRegion, double cutoffRadius,
			bool coversWholeDomain[3], const std::function<bool(const int[3])>& condition){
	std::vector<HaloRegion> regions;
		int d[3];
		for (d[0] = -1; d[0] <= 1; d[0]++) {
			for (d[1] = -1; d[1] <= 1; d[1]++) {
				for (d[2] = -1; d[2] <= 1; d[2]++) {
					if ((d[0] || d[1] || d[2]) == 0)  // if all are 0 (false), then continue
						continue;

					// we don't include anything, that does not conform with the condition
					if (!condition(d)) {
						continue;
					}

					HaloRegion tmp = initialRegion;
					for (unsigned int dimension = 0; dimension < 3; dimension++) {
						if (d[dimension] == 0) {
							tmp.rmin[dimension] = initialRegion.rmin[dimension];
							tmp.rmax[dimension] = initialRegion.rmax[dimension];
						} else if (d[dimension] == -1) { // LOWER
							tmp.rmin[dimension] = initialRegion.rmin[dimension];
							tmp.rmax[dimension] = initialRegion.rmin[dimension] + cutoffRadius;
						} else { //d[dimension==1 - UPPER
							tmp.rmin[dimension] = initialRegion.rmax[dimension] - cutoffRadius;
							tmp.rmax[dimension] = initialRegion.rmax[dimension];
						}
					}
					tmp.offset[0] = d[0];
					tmp.offset[1] = d[1];
					tmp.offset[2] = d[2];

					tmp.width = cutoffRadius;
					regions.push_back(tmp);
				}
			}
		}
		return regions;
}
