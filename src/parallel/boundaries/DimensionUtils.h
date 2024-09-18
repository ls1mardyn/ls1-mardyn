/*
 * DimensionUtils.h
 *
 *  Created on: 18 Sep 2024
 *      Author: amartyads
 */

#pragma once

#include <string>
#include "utils/Math.h"

namespace DimensionUtils {
/**
 * enum storing the axes and direction.
 *
 * The dimensions are POSX, NEGX, POSY, NEGY, POSZ and NEGZ.
 *
 * This is hardcoded for 3D, and ERROR is included for sanity checks.
 */
enum class DimensionType { POSX, NEGX, POSY, NEGY, POSZ, NEGZ, ERROR };

/* Check if a dimension is allowed. */
bool isDimensionNumericPermissible(int dim);

/* Convert a dimension from number to DimensionType, where x = +-1, y = +-2 and
 * z = +-3 */
DimensionType convertNumericToDimension(int dim);

/**
 * Convert an LS1 internal dimension representation int to DimensionType, where
 * x = 0, y = 1 and z = 2 Since this does not contain direction information, the
 * positive direction is returned.
 */
DimensionType convertLS1DimIndexToEnumPositive(int dim);

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
int convertEnumToLS1DimIndex(DimensionType dimension);


/* Returns the sign of a number, used for determining direction from a
 * dimension. */
inline int findSign(DimensionType dimension) {
  return ::findSign(convertDimensionToNumeric(dimension));
}

} // namespace DimensionUtils