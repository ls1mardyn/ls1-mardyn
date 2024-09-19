/*
 * BoundaryUtils.h
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#pragma once

#include <string>

/**
 * @brief Includes the BoundaryType enum and helper functions for boundary types.
 * @author Amartya Das Sharma
 */
namespace BoundaryUtils {

/**
 * @brief enum storing the types of global boundaries currently supported.
 *
 * The currently supported boundary types are
 * PERIODIC - periodic boundaries, using the default (communicate particle transfers with neighbours) behaviour
 * OUTFLOW - molecules exiting the boundary are deleted
 * REFLECTING - molecules exiting the boundary have their velocities reversed in the direction they are leaving (i.e.
 * reflected) ERROR - kept for sanity checks
 *
 * This can be extended if needed.
 */
enum class BoundaryType { PERIODIC, OUTFLOW, REFLECTING, ERROR };

/**
 * @brief Convert strings read from the XML input file into BoundaryType enum members.
 *
 * The conversion logic is as follows:
 * - Any string with "per" - BoundaryType::PERIODIC
 * - Any string with "ref" - BoundaryType::REFLECTING
 * - Any string with "out" - BoundaryType::OUTFLOW
 * An error occurs if none of the substrings are found.
 *
 * @param boundary string read from xml file
 * @return BoundaryType which is the corresponding BoundaryType enum member
 */
BoundaryType convertStringToBoundary(const std::string &boundary);

/**
 * @brief Converts BoundaryType enum members into strings.
 *
 * Used mainly for logging.
 *
 * @param boundary BoundaryType enum member
 * @return std::string corresponding to the enum member
 */
std::string convertBoundaryToString(BoundaryType boundary);

} // namespace BoundaryUtils
