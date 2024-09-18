/*
 * RegionUtils.h
 *
 *  Created on: 18 Sep 2024
 *      Author: amartyads
 */

#pragma once

#include <tuple>
#include "DimensionUtils.h"
#include "molecules/Molecule.h"

namespace RegionUtils {

/**
 * When given a domain delimited by givenRegionBegin and givenRegionEnd
 * and a DimensionType, and a requested regionWidth,  this function
 * returns a cropped box from the domain. This box is of size regionWidth
 * in the dimension specified, and the size of the domain otherwise.
 *
 * E.g. if the dimension given is +x, and the regionWidth is 3, this function
 * returns a box, which would contain every particle that is within 3
 * units from the +x boundary wall, in the original domain
 * demarcated by givenRegion{begin,end}.
 */
std::tuple<std::array<double, 3>, std::array<double, 3>>
getInnerRegionSlab(const std::array<double, 3> &givenRegionBegin,
               const std::array<double, 3> &givenRegionEnd,
               DimensionUtils::DimensionType dimension, double regionWidth);

/**
 * Checks if a given molecule is leaving the given region in the next timestep,
 * with added scalar velocity adjustment, but only in the given dimension.
 *
 * The molecule's position and velocity in the specified dimension is extracted
 * from the molecule itself, and the next position is calculated using the
 * velocity adjustment and the timestep length.
 *
 * When this function is called, it is usually expected that the forces have
 * been updated for the current timestep and the positions and velocities
 * aren't. Hence the nextStepVelAdjustment is provided to the caller function,
 * so that it can handle the change in velocity due to forces. This is not done
 * by the util function itself, to keep it as generic as possible.
 */
bool isMoleculeLeaving(const Molecule &molecule,
                       const std::array<double, 3> &regionBegin,
                       const std::array<double, 3> &regionEnd,
                       DimensionUtils::DimensionType dimension, double timestepLength,
                       double nextStepVelAdjustment);

/**
 * When given a domain delimited by givenRegionBegin and givenRegionEnd
 * and a DimensionType, and a requested regionWidth,  this function
 * returns a box outside the domain. This box is of size regionWidth
 * in the dimension specified, and the size of the domain plus regionwidth
 * otherwise.
 *
 * E.g. if the dimension given is +x, and the regionWidth is 3, this function
 * returns a box, which contains every ghost particle that is 3 units away
 * from +x.
 */
std::tuple<std::array<double, 3>, std::array<double, 3>>
getOuterRegionSlab(const std::array<double, 3> &givenRegionBegin,
               const std::array<double, 3> &givenRegionEnd,
               DimensionUtils::DimensionType dimension, 
               const std::array<double, 3> &regionWidth);

} // namespace RegionUtils