//
// Created by alex on 05.04.24.
//

#ifndef MARDYN_REGION_H
#define MARDYN_REGION_H

#include "utils/xmlfileUnits.h"

#include <array>

class Domain;

namespace Resolution{
	/**
 * Implemented AdResS runs on 3 different resolutions.
 * FullParticle -> molecules are represented with full detail
 * Hybrid -> molecules in transition area are FullParticle and CoarseGrain at the same time
 * CoarseGrain -> simplified molecule representation with lower detail
 * */
	enum ResolutionType{
		FullParticle = 0, Hybrid = 1, CoarseGrain = 2, ResolutionCount
	};

	/**
 * FPRegion is the Full Particle region. Currently it is a box starting in _low and goes until _high.
 * The hybrid area is a shell around the FP box. Each FPRegion can have different dimensions for the hybrid area.
 * */
	struct FPRegion {
		/**
		 * Defines a boundary region.
		 * H_FP is the region between the hybrid and full particle region.
		 * CG_H is the region between the hybrid and coarse grain region.
		 * */
		enum Intersection { H_FP = 0, CG_H = 1};

		//! @brief front left lower corner of FP
		std::array<double, 3> _low;
		//! @brief rear right upper corner of FP
		std::array<double, 3> _high;
		//! @brief front left lower corner of Hybrid
		std::array<double, 3> _lowHybrid;
		//! @brief rear right upper corner of Hybrid
		std::array<double, 3> _highHybrid;
		//! @brief center point of FP region
		std::array<double, 3> _center;
		//! @brief thickness of hybrid walls
		std::array<double, 3> _hybridDims;
		//! @brief dimensions of the FP region
		std::array<double, 3> _dim;

		/**
		 * @brief Constructs a FPRegion
		 * @param low is the left lower corner of the FP box
		 * @param high is the right higher corner of the FP box
		 * @param hDim is the width of the hybrid area shell around the FP box
		 * */
		FPRegion(const std::array<double,3>& low = {0,0,0}, const std::array<double, 3>& high = {0,0,0}, const std::array<double,3>& dims = {0,0,0})
				: _low(low), _high(high), _hybridDims(dims), _dim({0,0,0}) {}

		//! @brief gets the lower corner of hybrid region
		[[nodiscard]] const std::array<double, 3>& getLowHybrid() const { return _lowHybrid; }

		//! @brief gets the upper corner of hybrid region
		[[nodiscard]] const std::array<double, 3>& getHighHybrid() const { return _highHybrid; }

		//! @brief computes all fields according to input data low,high,dims
		void init();

		//! @brief initialized this FPRegion instance from XML data
		void readXML(XMLfileUnits &xmlconfig);

		/**
		 * @brief Computes the intersection point of the line that is created by connection
		 * the center of the region and the provided point with the either the inner or outer box of this region
		 * depending on the intersection type.
		 * @param point outside or inside of the box that is used for the computation
		 * @param inter Intersection::H_FP inner box, Intersection::CG_H outer box
		 * */
		[[nodiscard]] std::array<double, 3> computeIntersection(std::array<double,3> point, Intersection inter) const;

		/**
		 * @brief checks if the point is in the provided box defined by low and high
		 * @param point is the point to be checked
		 * @param low is the inclusive lower left front corner of the box
		 * @param high is the exclusive upper right rear corner of the box
		 * @returns true when point is in the box
		 * */
		static bool isInnerPoint(std::array<double,3> point, std::array<double, 3> low, std::array<double, 3> high);

		/**
		 * @brief checks if the given point is inside the specified region, which is either the FP or H box defined by this.
		 * When this region goes beyond the domain bounds, it should wrap around as the simulation uses periodic bounds.
		 * This method also checks for that by using the passed domain.
		 * @param domain to check periodic bounds
		 * @param region either FullParticle or Hybrid
		 * @param point to check
		 * */
		bool isInnerPointDomain(Domain* domain, ResolutionType region, std::array<double, 3> point) const;

		/**
		 * @brief checks if this FPRegion (the FullParticle area) is in the box defined by low and high.
		 * @param low lower corner of box
		 * @param high upper corner of box
		 * */
		[[maybe_unused]] bool isRegionInBox(std::array<double, 3> low, std::array<double, 3> high);
	};
};

#endif //MARDYN_REGION_H
