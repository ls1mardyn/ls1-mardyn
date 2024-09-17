/*
 * BoundaryHandler.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include <algorithm>

#include "BoundaryHandler.h"

#include "Simulation.h"
#include "integrators/Integrator.h"

#include "utils/Logger.h"

BoundaryUtils::BoundaryType BoundaryHandler::getGlobalWallType(
    BoundaryUtils::DimensionType dimension) const {
  return _boundaries.at(dimension);
}

BoundaryUtils::BoundaryType
BoundaryHandler::getGlobalWallType(int dimension) const {
  return getGlobalWallType(
      BoundaryUtils::convertLS1DimsToDimensionPos(dimension));
}

void BoundaryHandler::setGlobalWallType(BoundaryUtils::DimensionType dimension,
                                        BoundaryUtils::BoundaryType value) {
  if (dimension != BoundaryUtils::DimensionType::ERROR)
    _boundaries[dimension] = value;
}

void BoundaryHandler::setGlobalRegion(const double *start, const double *end) {
  for (short int i = 0; i < 3; i++) {
    _globalRegionStart[i] = start[i];
    _globalRegionEnd[i] = end[i];
  }
}

void BoundaryHandler::setLocalRegion(const double *start, const double *end) {
  for (short int i = 0; i < 3; i++) {
    _localRegionStart[i] = start[i];
    _localRegionEnd[i] = end[i];
  }
}

void BoundaryHandler::updateGlobalWallLookupTable() {
  _isGlobalWall[BoundaryUtils::DimensionType::POSX] =
      BoundaryUtils::isNearRel(_localRegionEnd[0], _globalRegionEnd[0]);
  _isGlobalWall[BoundaryUtils::DimensionType::NEGX] =
      BoundaryUtils::isNearRel(_localRegionStart[0], _globalRegionStart[0]);
  _isGlobalWall[BoundaryUtils::DimensionType::POSY] =
      BoundaryUtils::isNearRel(_localRegionEnd[1], _globalRegionEnd[1]);
  _isGlobalWall[BoundaryUtils::DimensionType::NEGY] =
      BoundaryUtils::isNearRel(_localRegionStart[1], _globalRegionStart[1]);
  _isGlobalWall[BoundaryUtils::DimensionType::POSZ] =
      BoundaryUtils::isNearRel(_localRegionEnd[2], _globalRegionEnd[2]);
  _isGlobalWall[BoundaryUtils::DimensionType::NEGZ] =
      BoundaryUtils::isNearRel(_localRegionStart[2], _globalRegionStart[2]);
}

bool BoundaryHandler::hasInvalidBoundary() const {
  return std::any_of(_boundaries.begin(), _boundaries.end(),[](const auto &keyVal) {
    const auto [dim, boundaryType] = keyVal;
    return boundaryType == BoundaryUtils::BoundaryType::ERROR;
  });
}

bool BoundaryHandler::hasNonPeriodicBoundary() const {
  return std::any_of(_boundaries.begin(), _boundaries.end(),[](const auto &keyVal) {
    const auto [dim, boundaryType] = keyVal;
    return boundaryType == BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL;
  });
}

bool BoundaryHandler::isGlobalWall(
    BoundaryUtils::DimensionType dimension) const {
  return _isGlobalWall.at(dimension);
}

bool BoundaryHandler::isGlobalWall(int dimension) const {
  return isGlobalWall(BoundaryUtils::convertLS1DimsToDimensionPos(dimension));
}

void BoundaryHandler::processGlobalWallLeavingParticles(
    ParticleContainer *moleculeContainer, double timestepLength) {
  const auto cutoff = moleculeContainer->getCutoff();
  for (auto const [currentDim, currentWallIsGlobalWall] : _isGlobalWall) {
    if (!currentWallIsGlobalWall)
      continue;

    switch (getGlobalWallType(currentDim)) {
    case BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL:
      // default behaviour
      break;

    case BoundaryUtils::BoundaryType::OUTFLOW:
      [[fallthrough]];
    case BoundaryUtils::BoundaryType::REFLECTING: {
      // create region by using getInnerBuffer()
      const auto [curWallRegionBegin, curWallRegionEnd] =
          BoundaryUtils::getInnerBuffer(_localRegionStart, _localRegionEnd,
                                        currentDim, cutoff);
      // grab an iterator from the converted coords
      const auto particlesInRegion = moleculeContainer->regionIterator(
          curWallRegionBegin.data(), curWallRegionEnd.data(),
          ParticleIterator::ONLY_INNER_AND_BOUNDARY);

      // iterate through all molecules
      for (auto it = particlesInRegion; it.isValid(); ++it) {
        Molecule curMolecule = *it;

        // Calculate the change in velocity, which the leapfrog method will
        // apply in the next velocity update to the dimension of interest.
        const int currentDimInt =
            BoundaryUtils::convertDimensionToLS1Dims(currentDim);
        const double halfTimestep = .5 * timestepLength;
        const double halfTimestepByMass = halfTimestep / it->mass();
        const double force = it->F(currentDimInt);
        const double nextStepVelAdjustment = halfTimestepByMass * force;

        // check if the molecule would leave the bounds
        if (BoundaryUtils::isMoleculeLeaving(
                curMolecule, curWallRegionBegin, curWallRegionEnd, currentDim,
                timestepLength, nextStepVelAdjustment)) {
          if (getGlobalWallType(currentDim) ==
              BoundaryUtils::BoundaryType::REFLECTING) {
            double currentVel = it->v(currentDimInt);
            // change the velocity in the dimension of interest such that when
            // the leapfrog integrator adds nextStepVelAdjustment in the next
            // velocity update, the final result ends up being the intended,
            // reversed velocity: -(currentVel+nextStepVelAdjustment)
            it->setv(currentDimInt, -currentVel - nextStepVelAdjustment -
                                        nextStepVelAdjustment);
          } else { // outflow, delete the particle if it would leave
            moleculeContainer->deleteMolecule(it, false);
          }
        }
      }
      break;
    }
    default:
      Log::global_log->error()
          << "Boundary type error! Received type not allowed!" << std::endl;
      mardyn_exit(1);
    }
  }
}

void BoundaryHandler::removeNonPeriodicHalos(
    ParticleContainer *moleculeContainer) {
  // get halo lengths in each dimension
  double buffers[] = {moleculeContainer->get_halo_L(0),
                      moleculeContainer->get_halo_L(1),
                      moleculeContainer->get_halo_L(2)};
  for (auto const [currentDim, currentWallIsGlobalWall] : _isGlobalWall) {
    if (!currentWallIsGlobalWall)
      continue;

    switch (getGlobalWallType(currentDim)) {
    case BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL:
      // default behaviour
      break;

    case BoundaryUtils::BoundaryType::OUTFLOW:
      [[fallthrough]];
    case BoundaryUtils::BoundaryType::REFLECTING: {
      // create region by using getOuterBuffer()
      auto const [curWallRegionBegin, curWallRegionEnd] =
          BoundaryUtils::getOuterBuffer(_localRegionStart, _localRegionEnd,
                                        currentDim, buffers);

      // grab an iterator from the converted coords
      auto particlesInRegion = moleculeContainer->regionIterator(
          curWallRegionBegin.data(), curWallRegionEnd.data(), ParticleIterator::ALL_CELLS);
      for (auto it = particlesInRegion; it.isValid(); ++it) {
        // delete all halo particles
        moleculeContainer->deleteMolecule(it, false);
      }
      break;
    }
    default:
      Log::global_log->error()
          << "Boundary type error! Received type not allowed!" << std::endl;
      mardyn_exit(1);
    }
  }
#ifndef MARDYN_AUTOPAS
  moleculeContainer->updateBoundaryAndHaloMoleculeCaches();
#endif
}
