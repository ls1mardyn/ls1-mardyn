/*
 * BoundaryHandler.h
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#pragma once

#include "BoundaryUtils.h"
#include "particleContainer/ParticleContainer.h"

#include <array>
#include <map>
#include <string>
#include <vector>

/**
 * Class to handle boundary conditions, namely leaving and halo particles.
 *
 * The objects of this class store the local and global bounds of the subdomain
 * in every process, and provide functions to deal with leaving particles, and
 * delete halo particles.
 *
 * The internal walls of the subdomain, touching other subdomains are 'local'
 * walls while the walls that are also the limits of the global domain are
 * 'global' walls.
 *
 * Since the behaviour of 'local' walls are unchanged, they are assigned to
 * 'PERIODIC'.
 */

class BoundaryHandler {
public:
  BoundaryHandler();

  /* Find the boundary type of a global wall for a particular dimension. */
  BoundaryUtils::BoundaryType getGlobalWallType(std::string dimension) const;

  /* Set the boundary type of a global wall for a particular dimension. */
  void setGlobalWallType(std::string dimension, BoundaryUtils::BoundaryType value);

  /* Find the boundary type of a global wall for a particular dimension. */
  BoundaryUtils::BoundaryType
  getGlobalWallType(BoundaryUtils::DimensionType dimension) const;

  /* Set the boundary type of a global wall for a particular dimension. */
  void setGlobalWallType(BoundaryUtils::DimensionType dimension,
                     BoundaryUtils::BoundaryType value);
  BoundaryUtils::BoundaryType getGlobalWallType(int dimension) const;

  /* Check if any of the global boundaries have invalid types. */
  bool hasInvalidBoundary() const;

  /**
   *  Check if any of the global boundaries are non-periodic.
   *
   *  This check helps bypass all boundary-related code if default behaviour
   * (all periodic boundaries) is expected.
   */
  bool hasNonPeriodicBoundary() const;

  /* Set bounds for global subdomain. */
  void setGlobalRegion(double *start, double *end);

  /* Set bounds for local subdomain. */
  void setLocalRegion(double *start, double *end);

  /* Set bounds for global subdomain. */
  void setGlobalRegion(std::array<double, 3> start, std::array<double, 3> end);

  /* Set bounds for local subdomain.*/
  void setLocalRegion(std::array<double, 3> start, std::array<double, 3> end);

  /**
   * Determine which walls in the local region are actually global walls.
   *
   * Should be called after changing global and local regions (typically after a
   * rebalance).
   */
  void findGlobalWallsInLocalRegion();

  /* Check if the local wall in a particular dimension is actually a global
   * wall. */
  bool isGlobalWall(BoundaryUtils::DimensionType dimension) const;

  /* Check if the local wall in a particular dimension is actually a global
   * wall. */
  bool isGlobalWall(int dimension) const;

  /**
   * Processes all particles that would leave the global domain.
   *
   * If a subdomain has no global walls, this function does nothing.
   * For every global wall, the function iterates through all particles that are
   * within one cutoff distance away from the wall. If these particles would leave the
   * global box in the next simulation, the following is done:
   *
   * PERIODIC - nop (default behaviour).
   * REFLECTING - The particle's velocity is reversed normal to the wall it's
   * leaving. 
   * OUTFLOW - The particle is deleted.
   */
  void processGlobalWallLeavingParticles(ParticleContainer *moleculeContainer,
                                         double timestepLength);

  /**
   * Processes all halo particles outside the global domain.
   *
   * If a subdomain has no global walls, this function does nothing.
   * For every global wall, the function iterates through all halo particles
   * that are within one cutoff distance away from the wall. The following is done for
   * each particle:
   *
   * PERIODIC - nop (default behaviour).
   * REFLECTING - The halo particle is deleted.
   * OUTFLOW - The halo particle is deleted.
   */
  void removeNonPeriodicHalos(ParticleContainer *moleculeContainer);

private:
  /* List of global boundary type per dimension. */
  std::map<BoundaryUtils::DimensionType, BoundaryUtils::BoundaryType>
      _boundaries;

  /* Lookup table to check if a local wall is also global. */
  std::map<BoundaryUtils::DimensionType, bool> _isGlobalWall;

  /* Global region start/end. */
  std::array<double, 3> _globalRegionStart, _globalRegionEnd;

  /* Local region start/end. */
  std::array<double, 3> _localRegionStart, _localRegionEnd;
};
