/*
 * BoundaryUtils.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include <algorithm>
#include <string>

#include "BoundaryUtils.h"

#include "Simulation.h"
#include "utils/Logger.h"

bool BoundaryUtils::isDimensionStringPermissible(std::string dimension) {
  return std::find(permissibleDimensionsString.begin(),
                   permissibleDimensionsString.end(),
                   dimension) != permissibleDimensionsString.end();
}

bool BoundaryUtils::isDimensionNumericPermissible(int dim) {
  return (dim >= -3 && dim <= 3 && dim != 0);
}

BoundaryUtils::DimensionType
BoundaryUtils::convertStringToDimension(std::string dimension) {
  if (!isDimensionStringPermissible(dimension)) {
    Log::global_log->error()
        << "Invalid dimension passed for enum conversion" << std::endl;
    mardyn_exit(1);
    return DimensionType::ERROR;
  }
  if (dimension == "+x")
    return DimensionType::POSX;
  if (dimension == "+y")
    return DimensionType::POSY;
  if (dimension == "+z")
    return DimensionType::POSZ;
  if (dimension == "-x")
    return DimensionType::NEGX;
  if (dimension == "-y")
    return DimensionType::NEGY;
  // if(dimension == "-z")
  return DimensionType::NEGZ;
}

BoundaryUtils::DimensionType BoundaryUtils::convertNumericToDimension(int dim) {
  if (!isDimensionNumericPermissible(dim)) {
    Log::global_log->error()
        << "Invalid dimension passed for enum conversion" << std::endl;
    mardyn_exit(1);
    return DimensionType::ERROR;
  }
  switch (findSign(dim)) {
  case -1:
    switch (dim) {
    case 1:
      return DimensionType::NEGX;
    case 2:
      return DimensionType::NEGY;
    default:
      return DimensionType::NEGZ; // case 3
    }
  case 1:
    switch (dim) {
    case 1:
      return DimensionType::POSX;
    case 2:
      return DimensionType::POSY;
    default:
      return DimensionType::POSZ; // case 3
    }
  default: // should never happen
    Log::global_log->error()
        << "Invalid dimension passed for enum conversion" << std::endl;
    mardyn_exit(1);
  }
  return DimensionType::ERROR;
}

BoundaryUtils::DimensionType
BoundaryUtils::convertLS1DimsToDimensionPos(int dim) {
  DimensionType toRet = DimensionType::ERROR;
  switch (dim) {
  case 0:
    toRet = DimensionType::POSX;
    break;
  case 1:
    toRet = DimensionType::POSY;
    break;
  default: // case 2:
    toRet = DimensionType::POSZ;
  }
  return toRet;
}

std::vector<BoundaryUtils::DimensionType>
BoundaryUtils::convertHaloOffsetToDimensionVector(int *offset) {
  std::vector<DimensionType> toRet;
  for (int i = 0; i < 3; i++) {
    if (offset[i] != 0) {
      int numeric = (i + 1) * offset[i];
      toRet.push_back(convertNumericToDimension(numeric));
    }
  }
  return toRet;
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
  int toReturn = convertDimensionToNumeric(dimension);
  return findSign(toReturn) * toReturn;
}

int BoundaryUtils::convertDimensionToLS1Dims(DimensionType dimension) {
  return convertDimensionToNumericAbs(dimension) - 1;
}

BoundaryUtils::BoundaryType
BoundaryUtils::convertStringToBoundary(std::string boundary) {
  if (boundary == "periodic")
    return BoundaryType::PERIODIC;
  if (boundary == "reflecting" || boundary == "reflective")
    return BoundaryType::REFLECTING;
  if (boundary == "outflow")
    return BoundaryType::OUTFLOW;
  Log::global_log->error()
      << "Invalid boundary type passed. Check your input file!" << std::endl;
  return BoundaryType::ERROR;
}

std::tuple<std::array<double, 3>, std::array<double, 3>>
BoundaryUtils::getInnerBuffer(const std::array<double, 3> givenRegionBegin,
                              const std::array<double, 3> givenRegionEnd,
                              DimensionType dimension, double regionWidth) {
  std::array<double, 3> returnRegionBegin, returnRegionEnd;
  for (int i = 0; i < 3; i++) {
    returnRegionBegin[i] = givenRegionBegin[i];
    returnRegionEnd[i] = givenRegionEnd[i];
  }

  int dimensionLS1 = convertDimensionToLS1Dims(dimension);

  switch (dimension) // can be done with findsign() too, but this is clearer
  {
  case DimensionType::POSX:
  case DimensionType::POSY:
  case DimensionType::POSZ:
    returnRegionBegin[dimensionLS1] =
        returnRegionEnd[dimensionLS1] - regionWidth;
    break;

  case DimensionType::NEGX:
  case DimensionType::NEGY:
  case DimensionType::NEGZ:
    returnRegionEnd[dimensionLS1] =
        returnRegionBegin[dimensionLS1] + regionWidth;
    break;

  default:
    Log::global_log->error()
        << "Invalid dimension passed for inner buffer calculation" << std::endl;
    mardyn_exit(1);
  }
  return std::make_tuple(returnRegionBegin, returnRegionEnd);
}

bool BoundaryUtils::isMoleculeLeaving(const Molecule molecule,
                                      const std::array<double, 3> regionBegin,
                                      const std::array<double, 3> regionEnd,
                                      DimensionType dimension,
                                      double timestepLength,
                                      double nextStepVelAdjustment) {
  int ls1dim = convertDimensionToLS1Dims(dimension);
  int direction = findSign(dimension);
  double newPos =
      molecule.r(ls1dim) +
      (timestepLength * (molecule.v(ls1dim) + nextStepVelAdjustment));
  if (newPos < regionBegin[ls1dim] && direction < 0)
    return true;
  if (newPos > regionEnd[ls1dim] && direction > 0)
    return true;
  return false;
}

std::tuple<std::array<double, 3>, std::array<double, 3>>
BoundaryUtils::getOuterBuffer(const std::array<double, 3> givenRegionBegin,
                              const std::array<double, 3> givenRegionEnd,
                              DimensionType dimension, double *regionWidth) {
  std::array<double, 3> returnRegionBegin, returnRegionEnd;
  for (int i = 0; i < 3; i++) {
    returnRegionBegin[i] = givenRegionBegin[i];
    returnRegionEnd[i] = givenRegionEnd[i];
  }

  int dimensionLS1 = convertDimensionToLS1Dims(dimension);

  int extraDim1 = dimensionLS1 == 0 ? 1 : 0;
  int extraDim2 = dimensionLS1 == 2 ? 1 : 2;

  returnRegionBegin[extraDim1] =
      returnRegionBegin[extraDim1] - regionWidth[extraDim1];
  returnRegionEnd[extraDim1] =
      returnRegionEnd[extraDim1] + regionWidth[extraDim1];

  returnRegionBegin[extraDim2] =
      returnRegionBegin[extraDim2] - regionWidth[extraDim2];
  returnRegionEnd[extraDim2] =
      returnRegionEnd[extraDim2] + regionWidth[extraDim2];

  switch (dimension) // can be done with findsign() too, but this is clearer
  {
  case DimensionType::POSX:
  case DimensionType::POSY:
  case DimensionType::POSZ:
    returnRegionBegin[dimensionLS1] = returnRegionEnd[dimensionLS1];
    returnRegionEnd[dimensionLS1] =
        returnRegionBegin[dimensionLS1] + regionWidth[dimensionLS1];
    break;

  case DimensionType::NEGX:
  case DimensionType::NEGY:
  case DimensionType::NEGZ:
    returnRegionEnd[dimensionLS1] = returnRegionBegin[dimensionLS1];
    returnRegionBegin[dimensionLS1] =
        returnRegionEnd[dimensionLS1] - regionWidth[dimensionLS1];
    break;

  default:
    Log::global_log->error()
        << "Invalid dimension passed for inner buffer calculation" << std::endl;
    mardyn_exit(1);
  }
  return std::make_tuple(returnRegionBegin, returnRegionEnd);
}
