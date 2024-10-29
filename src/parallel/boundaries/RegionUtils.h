/*
 * RegionUtils.h
 *
 *  Created on: 18 Sep 2024
 *      Author: amartyads
 */

#pragma once

#include "DimensionUtils.h"
#include "molecules/Molecule.h"
#include <tuple>

/**
 * @brief Includes helper functions for region-related computation during boundary handling.
 * @author Amartya Das Sharma
 *
 * Does not actually work on any molecule containers, hence can be used
 * in a general sense for other calculations.
 */
namespace RegionUtils {

/**
 * @brief When given a domain delimited by givenRegionBegin and givenRegionEnd
 * and a DimensionType, and a requested regionWidth,  this function
 * returns a cropped box from the domain. This box is of size regionWidth
 * in the dimension specified, and the size of the domain otherwise.
 *
 * E.g. if the dimension given is +x, and the regionWidth is 3, this function
 * returns a box, which would contain every particle that is within 3
 * units from the +x boundary wall, in the original domain
 * demarcated by givenRegion{begin,end}.
 *
 * @param givenRegionBegin the start of the region to be cropped
 * @param givenRegionEnd the end of the region to be cropped
 * @param dimension the dimension which is cropped i.e. the resulting slab will be regionWidth wide in this dimension,
 * and unchanged otherwise.
 * @param regionWidth the width of the slab
 * @return std::tuple<std::array<double, 3>, std::array<double, 3>> the corners of the resultant slab inside the given
 * region
 */
std::tuple<std::array<double, 3>, std::array<double, 3>> getInnerRegionSlab(
	const std::array<double, 3> &givenRegionBegin, const std::array<double, 3> &givenRegionEnd,
	DimensionUtils::DimensionType dimension, double regionWidth);

/**
 * @brief Checks if a given molecule is leaving the given region in the next timestep,
 * with added scalar velocity adjustment, but only in the given dimension.
 *
 * The molecule's position and velocity in the specified dimension is extracted
 * from the molecule itself, and the next position is calculated using the
 * velocity adjustment and the timestep length.
 *
 * If the molecule was already out of the box (and stay out of the box after position adjustment), this function returns
 * true. If the molecule somehow ends up entering the box after the position adjustment, the function returns false.
 *
 * When this function is called, it is usually expected that the forces have
 * been updated for the current timestep and the positions and velocities
 * aren't. Hence the nextStepVelAdjustment is provided to the caller function,
 * so that it can handle the change in velocity due to forces. This is not done
 * by the util function itself, to keep it as generic as possible.
 *
 * @param molecule the molecule to be checked
 * @param regionBegin the start of the region to be checked
 * @param regionEnd the end of the region to be checked
 * @param dimension the dimension in which the molecule's new position must be calculated, to check for leaving
 * @param timestepLength the length of the whole time step, after which the molecule's position is calculated
 * @param nextStepVelAdjustment a flat scalar velocity added to the velocity in the `dimension` component after
 * calculations
 * @return true if the molecule would leave the box after timestepLength
 * @return false if the molecule would not leave the box after timestepLength
 */
bool isMoleculeLeaving(const Molecule &molecule, const std::array<double, 3> &regionBegin,
					   const std::array<double, 3> &regionEnd, DimensionUtils::DimensionType dimension,
					   double timestepLength, double nextStepVelAdjustment);

/**
 * @brief When given a domain delimited by givenRegionBegin and givenRegionEnd
 * and a DimensionType, and a requested regionWidth,  this function
 * returns a box outside the domain. This box is of size regionWidth
 * in the dimension specified, and the size of the domain plus regionwidth
 * otherwise.
 *
 * E.g. if the dimension given is +x, and the regionWidth is 3, this function
 * returns a box, which contains every ghost particle that is 3 units away
 * from +x.
 *
 * @param givenRegionBegin the start of the region under consideration
 * @param givenRegionEnd the end of the region under consideration
 * @param dimension the dimension in which the box should be sharing a wall with the givenRegion
 * @param regionWidth the width of the slab
 * @return std::tuple<std::array<double, 3>, std::array<double, 3>> the corners of the resultant region outside the
 * given region
 */
std::tuple<std::array<double, 3>, std::array<double, 3>> getOuterRegionSlab(
	const std::array<double, 3> &givenRegionBegin, const std::array<double, 3> &givenRegionEnd,
	DimensionUtils::DimensionType dimension, const std::array<double, 3> &regionWidth);

} // namespace RegionUtils