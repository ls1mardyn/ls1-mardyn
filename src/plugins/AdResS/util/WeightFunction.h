//
// Created by alex on 05.04.24.
//

#ifndef MARDYN_WEIGHTFUNCTION_H
#define MARDYN_WEIGHTFUNCTION_H

#include "Region.h"
#include "particleContainer/adapter/vectorization/SIMD_DEFINITIONS.h"

#include <array>

namespace Weight {

	using function_t = double (*) (std::array<double,3>, const Resolution::FPRegion&);

	/**
     * @brief Weighting function for AdResS force computation.
     * Implementation computes axis intersection points and uses euclidean distance to determine the period of
     * the underlying cosine function.
     * @param r position of the site
     * @param region region of FP
     * */
	double euclid(std::array<double, 3> r, const Resolution::FPRegion& region);

	/**
	 * @brief Weighting function for AdResS force computation.
	 * Implementation computes axis intersection points and uses manhattan distance to determine the period of
	 * the underlying cosine function.
	 * @param r position of the site
	 * @param region region of FP
	 * */
	double manhattan(std::array<double, 3> r, const Resolution::FPRegion& region);

	/**
	 * @brief Weighting function for AdResS force computation.
	 * Implementation computes weight percentage for each component by checking where each component is in the hybrid region.
	 * All component weights are multiplied together.
	 * This results in a weight function, where each component does its own contribution.
	 * @param r position of the site
	 * @param region region of FP
	 * */
	double component(std::array<double, 3> r, const Resolution::FPRegion& region);

	/**
	 * @brief Weighting function for AdResS force computation.
	 * Implementation conceptually find the nearest point on the surface of the inner region in respect to r and
	 * computes the distance between r and this point.
	 * By doing so, this weight function treats the FPRegion as a box with rounded corners and edges.
	 * @param r position of the site
	 * @param region region of FP
	 * */
	double nearest(std::array<double, 3> r, const Resolution::FPRegion& region);

	/**
	 * @brief Weighting function for AdResS force computation.
	 * Implementation disables Hybrid region. Only for testing purposes, do not use in production.
	 * @param r position of the site
	 * @param region region of FP
	 * */
	double flat(std::array<double, 3> r, const Resolution::FPRegion& region);

	//! @brief Array containing function pointers to all weight function implementations
	constexpr std::array<function_t, 5> functions {euclid, manhattan, component, nearest, flat};

	vcp_inline MaskCalcVec is_inner_vec(const RealCalcVec &r0, const RealCalcVec &r1, const RealCalcVec &r2,
											   const RealCalcVec &region_low_0, const RealCalcVec &region_low_1,
											   const RealCalcVec &region_low_2, const RealCalcVec &region_high_0,
											   const RealCalcVec &region_high_1, const RealCalcVec &region_high_2) {
		const MaskCalcVec dim0 = (r0 < region_high_0) and (region_low_0 < r0);
		const MaskCalcVec dim1 = (r1 < region_high_1) and (region_low_1 < r1);
		const MaskCalcVec dim2 = (r2 < region_high_2) and (region_low_2 < r2);
		return dim0 and dim1 and dim2;
	}

	vcp_inline RealCalcVec nearest_vec(const RealCalcVec& r0, const RealCalcVec& r1, const RealCalcVec& r2,
									   const RealCalcVec& region_low_0, const RealCalcVec& region_low_1, const RealCalcVec& region_low_2,
									   const RealCalcVec& region_high_0, const RealCalcVec& region_high_1, const RealCalcVec& region_high_2,
									   const RealCalcVec& hybrid_dim_0, const RealCalcVec& hybrid_dim_1, const RealCalcVec& hybrid_dim_2,
									   const RealCalcVec &hybrid_low_0, const RealCalcVec &hybrid_low_1, const RealCalcVec &hybrid_low_2,
									   const RealCalcVec &hybrid_high_0, const RealCalcVec &hybrid_high_1, const RealCalcVec &hybrid_high_2) {
		static const auto zero = RealCalcVec::zero();
		const MaskCalcVec in_FP = is_inner_vec(r0, r1, r2, region_low_0, region_low_1, region_low_2, region_high_0, region_high_1, region_high_2);
		const MaskCalcVec in_H = is_inner_vec(r0, r1, r2, hybrid_low_0, hybrid_low_1, hybrid_low_2, hybrid_high_0, hybrid_high_1, hybrid_high_2);

		const RealCalcVec dd0 = RealCalcVec::max(RealCalcVec::max(region_low_0 - r0, r0 - region_high_0), RealCalcVec::zero());
		const RealCalcVec dd1 = RealCalcVec::max(RealCalcVec::max(region_low_1 - r1, r1 - region_high_1), RealCalcVec::zero());
		const RealCalcVec dd2 = RealCalcVec::max(RealCalcVec::max(region_low_2 - r2, r2 - region_high_2), RealCalcVec::zero());
		const RealCalcVec dist = RealCalcVec::sqrt(dd0*dd0 + dd1*dd1 + dd2*dd2);
		auto b0 = RealCalcVec::cvt_MaskVec_to_RealCalcVec(dd0 != zero);
		auto b1 = RealCalcVec::cvt_MaskVec_to_RealCalcVec(dd1 != zero);
		auto b2 = RealCalcVec::cvt_MaskVec_to_RealCalcVec(dd2 != zero);
		const RealCalcVec hDim = RealCalcVec::sqrt(
				b0 * b0 * (hybrid_dim_0 * hybrid_dim_0) +
				b1 * b1 * (hybrid_dim_1 * hybrid_dim_1) +
				b2 * b2 * (hybrid_dim_2 * hybrid_dim_2)
		);

		const MaskCalcVec outside_rounded = hDim < dist;

		// instead of cos we use a poly of deg 3 here as it is a quite good approx for cos in range 0 to pi/2
		// it does not matter too much what we pick, it just needs to be differentiable and have zero grad at the borders
		// f(x) = 2x³-3x²+1
		static const RealCalcVec n_one = RealCalcVec::set1(-1);
		static const RealCalcVec one = RealCalcVec::set1(1);
		static const RealCalcVec two = RealCalcVec::set1(2);
		static const RealCalcVec three = RealCalcVec::set1(3);

		const RealCalcVec x = n_one * RealCalcVec::max(n_one * dist / hDim, n_one);
		const RealCalcVec x2 = x * x;
		const RealCalcVec weight = two * x * x2 - three * x2 + one;

		// remove ones in H but outside rounded
		const auto real_out_round = RealCalcVec::cvt_MaskVec_to_RealCalcVec(outside_rounded);
		const RealCalcVec weight_rounded = weight * (one - real_out_round * real_out_round);
		// remove all outside of H
		const auto real_in_h = RealCalcVec::cvt_MaskVec_to_RealCalcVec(in_H);
		const RealCalcVec weight_H = weight_rounded * real_in_h * real_in_h;
		// set FP to 1
		const auto real_in_fp = RealCalcVec::cvt_MaskVec_to_RealCalcVec(in_FP);
		const auto real_in_fp_pos = real_in_fp * real_in_fp;
		const RealCalcVec weight_FP = weight_H * (one - real_in_fp_pos) + one * real_in_fp_pos;
				//weight_H * (real_in_fp_pos / weight + (one - real_in_fp_pos));
		return weight_FP;
	}
}

#endif //MARDYN_WEIGHTFUNCTION_H
