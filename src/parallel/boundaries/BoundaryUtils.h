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
 * PERIODIC_OR_LOCAL - periodic boundaries or local boundaries: both of which
 * use the default (communicate particle transfers with neighbours) behaviour
 * OUTFLOW - molecules exiting the boundary are deleted
 * REFLECTING - molecules exiting the boundary have their velocities
 * reversed in the direction they are leaving
 * ERROR - kept for sanity checks
 *
 * This can be extended if needed.
 */
enum class BoundaryType { PERIODIC_OR_LOCAL, OUTFLOW, REFLECTING, ERROR };

/**
 * enum storing the axes and direction.
 *
 * The dimensions are POSX, NEGX, POSY, NEGY, POSZ and NEGZ.
 *
 * This is hardcoded for 3D, and ERROR is included for sanity checks.
 */
enum class DimensionType { POSX, NEGX, POSY, NEGY, POSZ, NEGZ, ERROR };

/* List of all allowed dimensions in string format. */
const std::array<std::string, 6> permissibleDimensionsString = {
    "+x", "+y", "+z", "-x", "-y", "-z"};

/* Check if a dimension is allowed. */
bool isDimensionStringPermissible(std::string dimension);

/* Check if a dimension is allowed. */
bool isDimensionNumericPermissible(int dim);

/* Convert a string from permissibleDimensionsString into DimensionType */
DimensionType convertStringToDimension(std::string dimension);

/* Convert a dimension from number to DimensionType, where x = +-1, y = +-2 and
 * z = +-3 */
DimensionType convertNumericToDimension(int dim);

/**
 * Convert an LS1 internal dimension representation int to DimensionType, where
 * x = 0, y = 1 and z = 2 Since this does not contain direction information, the
 * positive direction is returned.
 */
DimensionType convertLS1DimsToDimensionPos(int dim);

/* Convert a DimensionType into a string from permissibleDimensionsString */
std::string convertDimensionToString(DimensionType dimension);

/**
 * Convert a DimensionType into a string from permissibleDimensionsString, and
 * remove directions Ex: POSX and NEGX both return "x"
 */
std::string convertDimensionToStringAbs(DimensionType dimension);

/* Convert a dimension from DimensionType to number, where x = +-1, y = +-2 and
 * z = +-3 */
int convertDimensionToNumeric(DimensionType dimension);

/**
 * Convert a DimensionType into an int, and remove directions
 * +x/-x returns 1, +y/-y returns 2, +z/-z returns 3
 */
int convertDimensionToNumericAbs(DimensionType dimension);

/**
 * Convert a DimensionType to LS1 internal dimension representation int, where x
 * = 0, y = 1 and z = 2 Since this does not contain direction information, both
 * POSX and NEGX return 0, for example.
 */
int convertDimensionToLS1Dims(DimensionType dimension);

/* Used to convert string read from the XML input file. */
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
bool isMoleculeLeaving(const Molecule &molecule,
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

/* Returns the sign of a number, used for determining direction from a
 * dimension. */
inline int findSign(int n) { return n < 0 ? -1 : 1; }

/* Returns the sign of a number, used for determining direction from a
 * dimension. */
inline int findSign(DimensionType dimension) {
  return findSign(convertDimensionToNumeric(dimension));
}

/**
 * Used for equality checks when comparing floats.
 *
 * Taken from AutoPas - src/autopas/utils/Math.h
 */
bool isNearRel(double a, double b, double maxRelativeDifference = 1e-9);

} // namespace BoundaryUtils
