//
// Created by alex on 05.04.24.
//

#ifndef MARDYN_WEIGHTFUNCTION_H
#define MARDYN_WEIGHTFUNCTION_H

#include "Region.h"
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
}

#endif //MARDYN_WEIGHTFUNCTION_H
