//
// Created by alex on 05.04.24.
//
#include "WeightFunction.h"
#include "Region.h"
#include <cmath>


double Weight::euclid(std::array<double, 3> r, const Resolution::FPRegion &region) {
	// if point is in the FP region -> weight is 1
	if(Resolution::FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
		// point is in hybrid region -> weight is between 0 and 1
	else if(Resolution::FPRegion::isInnerPoint(r, region._lowHybrid, region._highHybrid)) {
		std::array<double, 3> intersect_inner = region.computeIntersection(r, Resolution::FPRegion::Intersection::H_FP);
		std::array<double, 3> intersect_outer = region.computeIntersection(r, Resolution::FPRegion::Intersection::CG_H);
		double hyb_axis_length = sqrt(std::pow(intersect_outer[0]-intersect_inner[0], 2) +
									  std::pow(intersect_outer[1]-intersect_inner[1], 2) +
									  std::pow(intersect_outer[2]-intersect_inner[2], 2));
		double dist = sqrt(std::pow(r[0]-intersect_inner[0], 2) +
						   std::pow(r[1]-intersect_inner[1], 2) +
						   std::pow(r[2]-intersect_inner[2], 2));

		return std::pow(std::cos(M_PI/(2*hyb_axis_length) * dist), 2);
	}
		// point is in the CG region -> weight is 0
	else return 0.;
}

double Weight::manhattan(std::array<double, 3> r, const Resolution::FPRegion &region) {
	// if point is in the FP region -> weight is 1
	if(Resolution::FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
		// point is in hybrid region -> weight is between 0 and 1
	else if(Resolution::FPRegion::isInnerPoint(r, region._lowHybrid, region._highHybrid)) {
		std::array<double, 3> intersect_inner = region.computeIntersection(r, Resolution::FPRegion::Intersection::H_FP);
		std::array<double, 3> intersect_outer = region.computeIntersection(r, Resolution::FPRegion::Intersection::CG_H);
		double hyb_axis_length = intersect_outer[0]-intersect_inner[0]  +
								 intersect_outer[1]-intersect_inner[1]  +
								 intersect_outer[2]-intersect_inner[2];
		double dist = r[0]-intersect_inner[0] +
					  r[1]-intersect_inner[1] +
					  r[2]-intersect_inner[2];

		return std::pow(std::cos(M_PI/(2*hyb_axis_length) * dist), 2);
	}
		// point is in the CG region -> weight is 0
	else return 0.;
}

double Weight::component(std::array<double, 3> r, const Resolution::FPRegion &region) {
	// if point is in the FP region -> weight is 1
	if(Resolution::FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
		// point is in hybrid region -> weight is between 0 and 1
	else if(Resolution::FPRegion::isInnerPoint(r, region._lowHybrid, region._highHybrid)) {
		std::array<bool, 3> contribute_dir = {r[0] <= region._low[0] || r[0] >= region._high[0],
											  r[1] <= region._low[1] || r[1] >= region._high[1],
											  r[2] <= region._low[2] || r[2] >= region._high[2]};
		std::array<double,3> dist_dir{std::max(std::max(region._low[0] - r[0], 0.), std::max(r[0] - region._high[0], 0.)),
									  std::max(std::max(region._low[1] - r[1], 0.), std::max(r[1] - region._high[1], 0.)),
									  std::max(std::max(region._low[2] - r[2], 0.), std::max(r[2] - region._high[2], 0.))};
		double w = 1.;
		for(int d = 0; d < 3; d++) {
			if(region._hybridDims[d] != 0)
				w *= contribute_dir[d] * std::cos(M_PI/(2*region._hybridDims[d]) * dist_dir[d]) + (1 - contribute_dir[d]);
		}
		return w;
	}
		// point is in the CG region -> weight is 0
	else return 0.;
}

double Weight::nearest(std::array<double, 3> r, const Resolution::FPRegion &region) {
	// if point is in the FP region -> weight is 1
	if(Resolution::FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
		// point is in hybrid region -> weight is between 0 and 1
	else if(Resolution::FPRegion::isInnerPoint(r, region._lowHybrid, region._highHybrid)) {
		//example: dist_dir = {dx, 0, 0} if point is either left or right of region
		//         dist_dir = {dx, dy, 0} if point is in edged that go front to rear in block
		//         dist_dir = {dx, dy, dz} if point is in one corner
		std::array<double,3> dist_dir{std::max(std::max(region._low[0] - r[0], 0.), std::max(r[0] - region._high[0], 0.)),
									  std::max(std::max(region._low[1] - r[1], 0.), std::max(r[1] - region._high[1], 0.)),
									  std::max(std::max(region._low[2] - r[2], 0.), std::max(r[2] - region._high[2], 0.))};
		// => distance to the closest point of r on the surface of the inner region is the l2-norm of dist_dir
		double dist = sqrt(std::pow(dist_dir[0],2)+std::pow(dist_dir[1],2)+std::pow(dist_dir[2],2));
		// compute hybrid size as hybrid width is not equal on all sides
		double hDim = 0;
		for(int d = 0; d < 3; d++) {
			hDim += (dist_dir[d] != 0) * std::pow(region._hybridDims[d], 2);
		}
		hDim = sqrt(hDim);
		if(dist >= hDim) return 0.; // we are outside the rounded region -> treat as CG
		return std::pow(std::cos(M_PI/(2*hDim) * dist), 2);
	}
		// point is in the CG region -> weight is 0
	else return 0.;
}

double Weight::flat(std::array<double, 3> r, const Resolution::FPRegion &region) {
	if(Resolution::FPRegion::isInnerPoint(r, region._low, region._high)) return 1.;
	return 0.;
}
