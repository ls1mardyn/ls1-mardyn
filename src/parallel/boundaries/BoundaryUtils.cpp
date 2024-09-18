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

#include "Simulation.h"
#include "utils/Logger.h"

bool BoundaryUtils::isDimensionNumericPermissible(int dim) {
  // permissible dimensions are {-1, -2, -3, 1, 2, 3}
  return (dim >= -3 && dim <= 3 && dim != 0);
}

BoundaryUtils::DimensionType BoundaryUtils::convertNumericToDimension(int dim) {
  if (!isDimensionNumericPermissible(dim)) {
    Log::global_log->error()
        << "Invalid dimension passed for enum conversion. Received value: " << dim << std::endl;
    mardyn_exit(1);
    return DimensionType::ERROR;
  }
  switch (dim) {
    case 1:
      return DimensionType::POSX;
    case 2:
      return DimensionType::POSY;
    case 3:
      return DimensionType::POSZ; // case 3
    case -1:
      return DimensionType::NEGX;
    case -2:
        return DimensionType::NEGY;
    case -3:
        return DimensionType::NEGZ;
  }
  return DimensionType::ERROR; //warning suppression
}

BoundaryUtils::DimensionType
BoundaryUtils::convertLS1DimIndexToEnumPos(int dim) {
  switch (dim) {
  case 0:
    return DimensionType::POSX;
  case 1:
    return DimensionType::POSY;
  case 2:
    return DimensionType::POSZ;
  default:
    return DimensionType::ERROR;
  }
}

std::string BoundaryUtils::convertDimensionToString(DimensionType dimension) {
  switch (dimension) {
  case DimensionType::POSX:
    return "+x";

  case DimensionType::POSY:
    return "+y";

  case DimensionType::POSZ:
    return "+z";

  case DimensionType::NEGX:
    return "-x";

  case DimensionType::NEGY:
    return "-y";

  case DimensionType::NEGZ:
    return "-z";

  default:
    Log::global_log->error()
        << "Invalid dimension passed for enum conversion" << std::endl;
    mardyn_exit(1);
    return "error";
  }
}

std::string
BoundaryUtils::convertDimensionToStringAbs(DimensionType dimension) {
  return convertDimensionToString(dimension).substr(1, 1);
}

int BoundaryUtils::convertDimensionToNumeric(DimensionType dimension) {
  switch (dimension) {
    case DimensionType::POSX:
      return 1;
    case DimensionType::POSY:
      return 2;
    case DimensionType::POSZ:
      return 3;
    case DimensionType::NEGX:
      return -1;
    case DimensionType::NEGY:
      return -2;
    case DimensionType::NEGZ:
      return -3;
    default:
      Log::global_log->error()
          << "Invalid dimension passed for enum conversion" << std::endl;
      mardyn_exit(1);
      return 0;
  }
}

int BoundaryUtils::convertDimensionToNumericAbs(DimensionType dimension) {
  return std::abs(convertDimensionToNumeric(dimension));
}

int BoundaryUtils::convertEnumToLS1DimIndex(DimensionType dimension) {
  return convertDimensionToNumericAbs(dimension) - 1;
}

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

std::tuple<std::array<double, 3>, std::array<double, 3>>
BoundaryUtils::getInnerRegionSlab(const std::array<double, 3> &givenRegionBegin,
                              const std::array<double, 3> &givenRegionEnd,
                              DimensionType dimension, double regionWidth) {

  std::array<double, 3> returnRegionBegin = givenRegionBegin;
  std::array<double, 3> returnRegionEnd = givenRegionBegin;

  const int dimensionLS1 = convertEnumToLS1DimIndex(dimension);
  switch (dimension) {
  // in positive case, set the beginning to end-width, or whole domain if width
  // too large
  case DimensionType::POSX:
    [[fallthrough]];
  case DimensionType::POSY:
    [[fallthrough]];
  case DimensionType::POSZ:
    returnRegionBegin[dimensionLS1] =
        std::max(returnRegionEnd[dimensionLS1] - regionWidth,
                 givenRegionBegin[dimensionLS1]);
    break;
  // in negative case, set the end to beginning+width, or whole domain if width
  // too large
  case DimensionType::NEGX:
    [[fallthrough]];
  case DimensionType::NEGY:
    [[fallthrough]];
  case DimensionType::NEGZ:
    returnRegionEnd[dimensionLS1] =
        std::min(returnRegionBegin[dimensionLS1] + regionWidth,
                 givenRegionEnd[dimensionLS1]);
    break;

  default:
    Log::global_log->error()
        << "Invalid dimension passed for inner buffer calculation\n";
    mardyn_exit(1);
  }
  return {returnRegionBegin, returnRegionEnd};
}

bool BoundaryUtils::isMoleculeLeaving(const Molecule &molecule,
                                      const std::array<double, 3> &regionBegin,
                                      const std::array<double, 3> &regionEnd,
                                      DimensionType dimension,
                                      double timestepLength,
                                      double nextStepVelAdjustment) {
  const int ls1dim = convertEnumToLS1DimIndex(dimension);
  const int direction = findSign(dimension);
  const double newPos =
      molecule.r(ls1dim) +
      (timestepLength * (molecule.v(ls1dim) + nextStepVelAdjustment));
  if (newPos <= regionBegin[ls1dim] && direction < 0)
    return true;
  if (newPos >= regionEnd[ls1dim] && direction > 0)
    return true;
  return false;
}

std::tuple<std::array<double, 3>, std::array<double, 3>>
BoundaryUtils::getOuterRegionSlab(const std::array<double, 3> &givenRegionBegin,
                              const std::array<double, 3> &givenRegionEnd,
                              DimensionType dimension, 
                              const std::array<double, 3> &regionWidth) {
  std::array<double, 3> returnRegionBegin = givenRegionBegin;
  std::array<double, 3> returnRegionEnd = givenRegionEnd;

  for (int i = 0; i < 3; i++) {
    returnRegionBegin[i] = givenRegionBegin[i];
    returnRegionEnd[i] = givenRegionEnd[i];
  }

  const int dimensionLS1 = convertEnumToLS1DimIndex(dimension);

  // find the two dimensions that are not being considered
  const int extraDim1 = dimensionLS1 == 0 ? 1 : 0;
  const int extraDim2 = dimensionLS1 == 2 ? 1 : 2;

  // extend the extra dimensions to cover all ghost areas
  returnRegionBegin[extraDim1] -= regionWidth[extraDim1];
  returnRegionEnd[extraDim1] += regionWidth[extraDim1];

  returnRegionBegin[extraDim2] -= regionWidth[extraDim2];
  returnRegionEnd[extraDim2] += regionWidth[extraDim2];

  switch (dimension) // can be done with findsign() too, but this is clearer
  {
  // in positive case, move the box begin to edge of domain, and box end to
  // beyond
  case DimensionType::POSX:
    [[fallthrough]];
  case DimensionType::POSY:
    [[fallthrough]];
  case DimensionType::POSZ:
    returnRegionBegin[dimensionLS1] = returnRegionEnd[dimensionLS1];
    returnRegionEnd[dimensionLS1] += regionWidth[dimensionLS1];
    break;

  // in negative case, move the box end to edge of domain, and box begin to
  // beyond
  case DimensionType::NEGX:
    [[fallthrough]];
  case DimensionType::NEGY:
    [[fallthrough]];
  case DimensionType::NEGZ:
    returnRegionEnd[dimensionLS1] = returnRegionBegin[dimensionLS1];
    returnRegionBegin[dimensionLS1] -= regionWidth[dimensionLS1];
    break;

  default:
    Log::global_log->error()
        << "Invalid dimension passed for inner buffer calculation" << std::endl;
    mardyn_exit(1);
  }
  return std::make_tuple(returnRegionBegin, returnRegionEnd);
}
