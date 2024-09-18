/*
 * BoundaryUtils.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include <algorithm>
#include <string>
#include <cctype>

#include "BoundaryUtils.h"

#include "utils/mardyn_assert.h"
#include "utils/Logger.h"

BoundaryUtils::BoundaryType
BoundaryUtils::convertStringToBoundary(const std::string &boundary) {
  std::string boundaryLowercase;
  std::transform(boundary.begin(), boundary.end(), boundaryLowercase.begin(),
    [](unsigned char c){ return std::tolower(c); });
  if (boundaryLowercase.find("per") != std::string::npos)
    return BoundaryType::PERIODIC_OR_LOCAL;
  if (boundaryLowercase.find("ref") != std::string::npos)
    return BoundaryType::REFLECTING;
  if (boundaryLowercase.find("out") != std::string::npos)
    return BoundaryType::OUTFLOW;
  Log::global_log->error()
      << "Invalid boundary type passed to "
         "BoundaryUtils::convertStringToBoundary. Check your input file!"
      << std::endl;
  mardyn_exit(1);
  return BoundaryType::ERROR; // warning suppression
}

std::string BoundaryUtils::convertBoundaryToString(BoundaryType boundary) {
  switch (boundary) {
  case BoundaryType::PERIODIC_OR_LOCAL:
    return "periodic";
  case BoundaryType::REFLECTING:
    return "reflecting";
  case BoundaryType::OUTFLOW:
    return "outflow";
  default:
    Log::global_log->error()
        << "Invalid boundary type passed to "
           "BoundaryUtils::convertBoundaryToString. Check your input file!"
        << std::endl;
    mardyn_exit(1);
  }
  return "error"; // warning suppression
}
