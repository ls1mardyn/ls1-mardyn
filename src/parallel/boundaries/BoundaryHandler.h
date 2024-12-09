/*
 * BoundaryHandler.h
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#pragma once

#include "BoundaryUtils.h"
#include "DimensionUtils.h"
#include "RegionUtils.h"
#include "particleContainer/ParticleContainer.h"

#include <array>
#include <map>
#include <string>
#include <vector>

/**
 * @brief Class to handle boundary conditions, namely leaving and halo particles.
 * @author Amartya Das Sharma
 *
 * The objects of this class store the local and global bounds of the subdomain in every process, and provide functions
 * to deal with leaving particles, and delete halo particles.
 *
 * The internal walls of the subdomain, touching other subdomains are 'local' walls while the walls that are also the
 * limits of the global domain are 'global' walls.
 *
 * All subdomains have a copy of the global wall types, and do not store their local boundary conditions; instead they
 * check whether a particular local wall is also a global wall, and then use the global lookup table to ascertain what
 * boundary effects to use.
 *
 * The default state for all global walls is 'PERIODIC', since that is the default ls1 behaviour.
 */
class BoundaryHandler {
public:
	BoundaryHandler() = default;

	/**
	 * @brief Find the boundary type of a global wall for a particular dimension.
	 *
	 * Since this returns the global wall type, every subdomain on every rank will return
	 * the same value for the same input.
	 *
	 * @param dimension the dimension of interest, must be DimensionType enum member
	 * @return BoundaryUtils::BoundaryType enum member corresponding to the global wall type
	 */
	BoundaryUtils::BoundaryType getGlobalWallType(DimensionUtils::DimensionType dimension) const;

	/**
	 * @brief Find the boundary type of a global wall for a particular dimension.
	 *
	 * Since this returns the global wall type, every subdomain on every rank will return
	 * the same value for the same input.
	 *
	 * @param dimension the dimension of interest, must be an ls1-style index (0 for x, 1 for y, 2 for z)
	 * @return BoundaryUtils::BoundaryType enum member corresponding to the global wall type
	 */
	BoundaryUtils::BoundaryType getGlobalWallType(int dimension) const;

	/**
	 * @brief Set the boundary type of a global wall for a particular dimension.
	 *
	 * @param dimension the dimension of interest, must be DimensionType enum member
	 * @param value the type of the boundary, must be BoundaryType enum member
	 */
	void setGlobalWallType(DimensionUtils::DimensionType dimension, BoundaryUtils::BoundaryType value);

	/**
	 * @brief Check if any of the global boundaries have invalid types.
	 *
	 * @return true if any of the global boundaries is BoundaryType::ERROR
	 * @return false if none of the global boundaries are BoundaryType::ERROR
	 */
	bool hasGlobalInvalidBoundary() const;

	/**
	 * @brief Check if any of the global boundaries are non-periodic.
	 *
	 * This check helps bypass all boundary-handling related code if default behaviour
	 * (all periodic boundaries) is found.
	 *
	 * @return true if any of the global boundaries are non-periodic
	 * @return false if all of the global boundaries are periodic
	 */
	bool hasGlobalNonPeriodicBoundary() const;

	/**
	 * @brief Set bounds for global domain.
	 *
	 * @param start double[3] array with coordinates of the start point
	 * @param end double[3] array with coordinates of the end point
	 */
	void setGlobalRegion(const double *start, const double *end);

	/**
	 * @brief Set bounds for local subdomain.
	 *
	 * @param start double[3] array with coordinates of the start point
	 * @param end double[3] array with coordinates of the end point
	 */
	void setLocalRegion(const double *start, const double *end);

	/**
	 * @brief Determine which walls in the local region are actually global walls.
	 *
	 * Should be called after changing global and local regions (typically after a rebalance).
	 */
	void updateGlobalWallLookupTable();

	/**
	 * @brief Check if the local wall in a particular dimension is actually a global wall.
	 *
	 * @param dimension the dimension of interest, must be DimensionType enum member
	 * @return true if the local wall in the given dimension is a global wall
	 * @return false if the local wall in the given dimension is not a global wall
	 */
	bool isGlobalWall(DimensionUtils::DimensionType dimension) const;

	/**
	 * @brief Check if the local wall in a particular dimension is actually a global wall.
	 *
	 * @param dimension the dimension of interest, must be an ls1-style index (0 for x, 1 for y, 2 for z)
	 * @return true if the local wall in the given dimension is a global wall
	 * @return false if the local wall in the given dimension is not a global wall
	 */
	bool isGlobalWall(int dimension) const;

	/**
	 * @brief Processes all particles that would leave the global domain.
	 *
	 * If a subdomain has no global walls, this function does nothing.
	 * For every global wall, the function iterates through all particles that are
	 * within one cutoff distance away from the wall. If these particles would
	 * leave the global box in the next simulation, the following is done:
	 *
	 * PERIODIC - No actions taken (default behaviour).
	 * REFLECTING - The particle's velocity is reversed normal to the wall it's leaving.
	 * OUTFLOW - The particle is deleted.
	 *
	 * If a particle would hit multiple walls with multiple types, the following happens:
	 * - If any of the walls is an outflow wall, the particle is deleted immediately, otherwise the following happen
	 * - The particle's velocity is reversed for the components where the wall is reflecting
	 * - Periodic effects happen as normal
	 *
	 * @param moleculeContainer used to get the cutoff, region iterators, and to delete particles if needed
	 * @param timestepLength the full timestep length, used to determine the position of the molecule
	 * in the next timestep to check whether it is leaving
	 */
	void processGlobalWallLeavingParticles(ParticleContainer *moleculeContainer, double timestepLength) const;

	/**
	 * @brief Processes all halo particles outside the global domain.
	 *
	 * If a subdomain has no global walls, this function does nothing.
	 * For every global wall, the function iterates through all halo particles
	 * that are within one cutoff distance away from the wall. The following is
	 * done for each particle:
	 *
	 * PERIODIC - No actions taken (default behaviour).
	 * REFLECTING / OUTFLOW - The halo particle is deleted, so that particles
	 * approaching the boundary do not decelerate due to influence from the halo
	 * particles, and preserve their velocities before being bounced/deleted
	 *
	 * @param moleculeContainer used to get the cutoff, region iterators, and to delete particles if needed
	 */
	void removeNonPeriodicHalos(ParticleContainer *moleculeContainer) const;

private:
	/**
	 * @brief Lookup table for global boundary type in every dimension.
	 *
	 * Set as periodic by default, since the default behaviour of LS1 is all periodic boundaries.
	 * This allows the <boundaries> tag to be optional.
	 */
	std::map<DimensionUtils::DimensionType, BoundaryUtils::BoundaryType> _boundaries{
		{DimensionUtils::DimensionType::POSX, BoundaryUtils::BoundaryType::PERIODIC},
		{DimensionUtils::DimensionType::POSY, BoundaryUtils::BoundaryType::PERIODIC},
		{DimensionUtils::DimensionType::POSZ, BoundaryUtils::BoundaryType::PERIODIC},
		{DimensionUtils::DimensionType::NEGX, BoundaryUtils::BoundaryType::PERIODIC},
		{DimensionUtils::DimensionType::NEGY, BoundaryUtils::BoundaryType::PERIODIC},
		{DimensionUtils::DimensionType::NEGZ, BoundaryUtils::BoundaryType::PERIODIC}};

	/* Lookup table to check if a local wall is also global. */
	std::map<DimensionUtils::DimensionType, bool> _isGlobalWall;

	/* Global region start/end. */
	std::array<double, 3> _globalRegionStart, _globalRegionEnd;

	/* Local region start/end. */
	std::array<double, 3> _localRegionStart, _localRegionEnd;
};
