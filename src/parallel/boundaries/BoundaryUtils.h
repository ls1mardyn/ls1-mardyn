/*
 * BoundaryUtils.h
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#pragma once

#include <array>
#include <string>
#include <tuple>
#include <vector>

#include "molecules/Molecule.h"

/**
 * Includes enums and helper functions for processing boundary conditions.
 *
 * Does not actually work on any molecule containers, hence can be used
 * in a general sense for other calculations.
 */
namespace BoundaryUtils {

/**
 * enum storing the types of boundaries currently supported.
 *
 * The currently supported boundary types are
 * PERIODIC - default behaviour, periodic boundaries
 * OUTFLOW - molecules exiting the boundary are deleted
 * REFLECTING - molecules exiting the boundary have their velocities
 * reversed in the direction they are leaving
 * ERROR - kept for sanity checks
 *
 * This can be extended if needed.
 */
enum class BoundaryType { PERIODIC, OUTFLOW, REFLECTING, ERROR };

/**
 * enum storing the axes and direction.
 *
 * This is fixed for 3D, and ERROR is included for sanity checks
 */
enum class DimensionType { POSX, NEGX, POSY, NEGY, POSZ, NEGZ, ERROR };

const std::array<std::string, 6> permissibleDimensionsString = {
    "+x", "+y", "+z", "-x", "-y", "-z"};
const std::array<int, 6> permissibleDimensionsInt = {-1, -2, -3, 1, 2, 3};

bool isDimensionStringPermissible(std::string dimension);
bool isDimensionNumericPermissible(int dim);

DimensionType convertStringToDimension(std::string dimension);
DimensionType convertNumericToDimension(int dim);
DimensionType convertLS1DimsToDimensionPos(int dim);
std::vector<DimensionType> convertHaloOffsetToDimensionVector(int *offset);

std::string convertDimensionToString(DimensionType dimension);
std::string convertDimensionToStringAbs(DimensionType dimension);
int convertDimensionToNumeric(DimensionType dimension);
int convertDimensionToNumericAbs(DimensionType dimension);
int convertDimensionToLS1Dims(DimensionType dimension);

/**
 * Used to convert string read from the XML input file
 */
BoundaryType convertStringToBoundary(std::string boundary);

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
getInnerBuffer(const std::array<double, 3> givenRegionBegin,
               const std::array<double, 3> givenRegionEnd,
               DimensionType dimension, double regionWidth);

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
bool isMoleculeLeaving(const Molecule molecule,
                       const std::array<double, 3> regionBegin,
                       const std::array<double, 3> regionEnd,
                       DimensionType dimension, double timestepLength,
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
getOuterBuffer(const std::array<double, 3> givenRegionBegin,
               const std::array<double, 3> givenRegionEnd,
               DimensionType dimension, double *regionWidth);

inline int findSign(int n) { return n < 0 ? -1 : 1; }
inline int findSign(DimensionType dimension) {
  return findSign(convertDimensionToNumeric(dimension));
}
} // namespace BoundaryUtils