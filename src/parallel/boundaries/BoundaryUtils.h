/*
 * BoundaryUtils.h
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#pragma once

#include <string>
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

/* Used to convert string read from the XML input file. */
BoundaryType convertStringToBoundary(const std::string &boundary);

std::string convertBoundaryToString(BoundaryType boundary);

} // namespace BoundaryUtils
