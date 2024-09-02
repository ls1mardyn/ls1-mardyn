/*
 * BoundaryHandler.h
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#pragma once

#include "BoundaryUtils.h"

#include <array>
#include <map>
#include <string>
#include <vector>

class BoundaryHandler {
public:
  BoundaryHandler();
  BoundaryUtils::BoundaryType getGlobalWall(std::string dimension) const;
  void setGlobalWall(std::string dimension, BoundaryUtils::BoundaryType value);
  BoundaryUtils::BoundaryType
  getGlobalWall(BoundaryUtils::DimensionType dimension) const;
  void setGlobalWall(BoundaryUtils::DimensionType dimension,
                     BoundaryUtils::BoundaryType value);
  BoundaryUtils::BoundaryType getGlobalWall(int dimension) const;
  bool hasInvalidBoundary() const;
  bool hasNonPeriodicBoundary() const;

  void setGlobalRegion(double *start, double *end);
  void setLocalRegion(double *start, double *end);

  void setGlobalRegion(std::array<double, 3> start, std::array<double, 3> end);
  void setLocalRegion(std::array<double, 3> start, std::array<double, 3> end);

  void findOuterWallsInLocalRegion();
  bool isOuterWall(BoundaryUtils::DimensionType dimension) const;
  bool isOuterWall(int dimension) const;
  bool processOuterWallLeavingParticles();
  void removeNonPeriodicHalos();

private:
  std::map<BoundaryUtils::DimensionType, BoundaryUtils::BoundaryType>
      _boundaries;
  std::map<BoundaryUtils::DimensionType, bool> _isOuterWall;
  std::array<double, 3> _globalRegionStart, _globalRegionEnd;
  std::array<double, 3> _localRegionStart, _localRegionEnd;
};