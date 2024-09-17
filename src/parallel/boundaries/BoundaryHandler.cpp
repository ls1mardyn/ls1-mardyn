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
#include "utils/Math.h"

#include "utils/Logger.h"

BoundaryUtils::BoundaryType BoundaryHandler::getGlobalWallType(
    BoundaryUtils::DimensionType dimension) const {
  return _boundaries.at(dimension);
}

BoundaryUtils::BoundaryType
BoundaryHandler::getGlobalWallType(int dimension) const {
  return getGlobalWallType(
      BoundaryUtils::convertLS1DimIndexToEnumPos(dimension));
}

void BoundaryHandler::setGlobalWallType(BoundaryUtils::DimensionType dimension,
                                        BoundaryUtils::BoundaryType value) {
  if (dimension != BoundaryUtils::DimensionType::ERROR)
    _boundaries[dimension] = value;
  else {
    Log::global_log->error()
          << "DimensionType::ERROR received in setGlobalWallType!!" << std::endl;
      mardyn_exit(1);
  }
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
      isNearRel(_localRegionEnd[0], _globalRegionEnd[0]);
  _isGlobalWall[BoundaryUtils::DimensionType::NEGX] =
      isNearRel(_localRegionStart[0], _globalRegionStart[0]);
  _isGlobalWall[BoundaryUtils::DimensionType::POSY] =
      isNearRel(_localRegionEnd[1], _globalRegionEnd[1]);
  _isGlobalWall[BoundaryUtils::DimensionType::NEGY] =
      isNearRel(_localRegionStart[1], _globalRegionStart[1]);
  _isGlobalWall[BoundaryUtils::DimensionType::POSZ] =
      isNearRel(_localRegionEnd[2], _globalRegionEnd[2]);
  _isGlobalWall[BoundaryUtils::DimensionType::NEGZ] =
      isNearRel(_localRegionStart[2], _globalRegionStart[2]);
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
  return isGlobalWall(BoundaryUtils::convertLS1DimIndexToEnumPos(dimension));
}

void BoundaryHandler::processGlobalWallLeavingParticles(
    ParticleContainer *moleculeContainer, double timestepLength) const {
  const auto cutoff = moleculeContainer->getCutoff();
  for (auto const [currentDim, currentWallIsGlobalWall] : _isGlobalWall) {
    if (!currentWallIsGlobalWall)
      continue;

    switch (getGlobalWallType(currentDim)) {
    case BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL:
      // nothing changes from normal ls1 behaviour, so leaving particles not touched by BoundaryHandler and are processed by
      // DomainDecompBase::handleDomainLeavingParticles()
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
      for (auto moleculeIter = particlesInRegion; moleculeIter.isValid(); ++moleculeIter) {
        // Calculate the change in velocity, which the leapfrog method will
        // apply in the next velocity update to the dimension of interest.
        const int currentDimInt =
            BoundaryUtils::convertEnumToLS1DimIndex(currentDim);
        const double halfTimestep = .5 * timestepLength;
        const double halfTimestepByMass = halfTimestep / moleculeIter->mass();
        const double force = moleculeIter->F(currentDimInt);
        const double nextStepVelAdjustment = halfTimestepByMass * force;

        // check if the molecule would leave the bounds
        if (BoundaryUtils::isMoleculeLeaving(
                *moleculeIter, curWallRegionBegin, curWallRegionEnd, currentDim,
                timestepLength, nextStepVelAdjustment)) {
          if (getGlobalWallType(currentDim) ==
              BoundaryUtils::BoundaryType::REFLECTING) {
            const double currentVel = moleculeIter->v(currentDimInt);
            // change the velocity in the dimension of interest such that when
            // the leapfrog integrator adds nextStepVelAdjustment in the next
            // velocity update, the final result ends up being the intended,
            // reversed velocity: -(currentVel+nextStepVelAdjustment)
            moleculeIter->setv(currentDimInt, -currentVel - nextStepVelAdjustment -
                                        nextStepVelAdjustment);
          } else { // outflow, delete the particle if it would leave
            moleculeContainer->deleteMolecule(moleculeIter, false);
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
    ParticleContainer *moleculeContainer) const {
  // get halo lengths in each dimension
  const std::array<double, 3> buffers = {
                      moleculeContainer->getHaloWidthForDimension(0),
                      moleculeContainer->getHaloWidthForDimension(1),
                      moleculeContainer->getHaloWidthForDimension(2)};
  for (auto const [currentDim, currentWallIsGlobalWall] : _isGlobalWall) {
    if (!currentWallIsGlobalWall)
      continue;

    switch (getGlobalWallType(currentDim)) {
    case BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL:
      // nothing changes from normal ls1 behaviour, so empty case, and halo particles left untouched
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
      for (auto moleculeIter = particlesInRegion; moleculeIter.isValid(); ++moleculeIter) {
        // delete all halo particles
        moleculeContainer->deleteMolecule(moleculeIter, false);
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
