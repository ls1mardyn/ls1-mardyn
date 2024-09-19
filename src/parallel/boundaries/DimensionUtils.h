/*
 * DimensionUtils.h
 *
 *  Created on: 18 Sep 2024
 *      Author: amartyads
 */

#pragma once

#include "utils/Math.h"
#include <string>

/**
 * @brief Includes the DimensionType enum and helper functions for dimension types.
 * @author Amartya Das Sharma
 */
namespace DimensionUtils {

/**
 * @brief enum storing the axes and direction.
 *
 * The dimensions are POSX, NEGX, POSY, NEGY, POSZ and NEGZ.
 *
 * This is hardcoded for 3D, and ERROR is included for sanity checks.
 */
enum class DimensionType { POSX, NEGX, POSY, NEGY, POSZ, NEGZ, ERROR };

/**
 * @brief Checks if a numeric dimension is allowed.
 *
 * Allowed numeric dimensions are +-1 for x, +-2 for y, +-3 for z
 *
 * @param dim integer dimension
 * @return true if `dim` is one of (1, 2, 3, -1, -2, -3)
 * @return false if `dim` is not one of (1, 2, 3, -1, -2, -3)
 */
bool isDimensionNumericPermissible(int dim);

/**
 * @brief Converts a dimension from number to DimensionType, where x = +-1, y = +-2 and z = +-3.
 *
 * Throws an error and exits if the dimension is not permissible.
 *
 * @param dim integer dimension
 * @return DimensionType enum member corresponding to the integer dimension
 */
DimensionType convertNumericToDimension(int dim);

/**
 * @brief Converts an LS1 internal dimension representation int to DimensionType, where
 * x = 0, y = 1 and z = 2 Since this does not contain direction information, the
 * positive direction is returned.
 *
 * Throws an error and exits if the dimension is not permissible.
 *
 * @param dim integer dimension (either 0, 1 or 2)
 * @return DimensionType enum member corresponding to the integer dimension
 */
DimensionType convertLS1DimIndexToEnumPositive(int dim);

/**
 * @brief Converts a DimensionType into a string.
 *
 * The expected return values are +x, -x, +y, -y, +z, -z. Can be used for logging purposes.
 * Throws an error and exits if DimensionType::ERROR is encountered.
 *
 * @param dimension DimensionType enum member
 * @return std::string which can be +x, -x, +y, -y, +z, -z
 */
std::string convertDimensionToString(DimensionType dimension);

/**
 * @brief Converts a DimensionType into a string, and remove directional information.
 *
 * POSX and NEGX both return "x", POSY and NEGY return "y", and POSZ and NEGZ return "z".
 * Can be used for logging purposes. Throws an error and exits if DimensionType::ERROR is encountered.
 *
 * @param dimension DimensionType enum member
 * @return std::string which can be x, y, or z
 */
std::string convertDimensionToStringAbs(DimensionType dimension);

/**
 * @brief Converts a dimension from DimensionType to number, where x = +-1, y = +-2 and
 * z = +-3. Throws an error and exits if DimensionType::ERROR is encountered.
 *
 * @param dimension DimensionType enum member
 * @return int which is one of (1, 2, 3, -1, -2, -3)
 */
int convertDimensionToNumeric(DimensionType dimension);

/**
 * @brief Converts a DimensionType into an int, and remove directions.
 *
 * +x/-x returns 1, +y/-y returns 2, +z/-z returns 3.
 * Throws an error and exits if DimensionType::ERROR is encountered.
 *
 * @param dimension DimensionType enum member
 * @return int which is one of (1, 2, 3)
 */
int convertDimensionToNumericAbs(DimensionType dimension);

/**
 * @brief Converts a DimensionType to LS1 internal dimension representation int, where x
 * = 0, y = 1 and z = 2
 *
 * Since this does not contain direction information, both POSX and NEGX return 0, for example.
 * Throws an error and exits if DimensionType::ERROR is encountered.
 *
 * @param dimension DimensionType enum member
 * @return int which is one of (0, 1, 2)
 */
int convertEnumToLS1DimIndex(DimensionType dimension);

/**
 * @brief Returns the direction from a given DimensionType.
 *
 * @param dimension DimensionType enum member
 * @return int which is 1 if the direction is positive, and -1 otherwise
 */
inline int findSign(DimensionType dimension) { return ::findSign(convertDimensionToNumeric(dimension)); }

} // namespace DimensionUtils