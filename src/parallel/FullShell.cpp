/*
 * FullShell.cpp
 *
 *  Created on: Oct 10, 2016
 *      Author: seckler
 */

#include "FullShell.h"

FullShell::FullShell() {
}

FullShell::~FullShell() {
}

std::vector<HaloRegion> FullShell::getHaloRegions(HaloRegion initialRegion, double cutoffRadius,
		bool coversWholeDomain[3]) {
	std::vector<HaloRegion> regions;
	int d[3];
	for (d[0] = -1; d[0] <= 1; d[0]++) {
		for (d[1] = -1; d[1] <= 1; d[1]++) {
			for (d[2] = -1; d[2] <= 1; d[2]++) {
				if ((d[0] == d[1]) == (d[2] == 0))
					continue;
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
				regions.push_back(tmp);
			}
		}
	}
	return regions;
}
