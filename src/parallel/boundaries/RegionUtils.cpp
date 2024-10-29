/*
 * RegionUtils.cpp
 *
 *  Created on: 18 Sep 2024
 *      Author: amartyads
 */

#include "RegionUtils.h"
#include "utils/mardyn_assert.h" //for mardyn_exit()

std::tuple<std::array<double, 3>, std::array<double, 3>> RegionUtils::getInnerRegionSlab(
	const std::array<double, 3> &givenRegionBegin, const std::array<double, 3> &givenRegionEnd,
	DimensionUtils::DimensionType dimension, double regionWidth) {

	std::array<double, 3> returnRegionBegin = givenRegionBegin;
	std::array<double, 3> returnRegionEnd = givenRegionEnd;

	const int dimensionLS1 = convertEnumToLS1DimIndex(dimension);
	switch (dimension) {
		// in positive case, set the beginning to end-width, or whole domain if width
		// too large
		case DimensionUtils::DimensionType::POSX:
			[[fallthrough]];
		case DimensionUtils::DimensionType::POSY:
			[[fallthrough]];
		case DimensionUtils::DimensionType::POSZ:
			returnRegionBegin[dimensionLS1] =
				std::max(returnRegionEnd[dimensionLS1] - regionWidth, givenRegionBegin[dimensionLS1]);
			break;
		// in negative case, set the end to beginning+width, or whole domain if width
		// too large
		case DimensionUtils::DimensionType::NEGX:
			[[fallthrough]];
		case DimensionUtils::DimensionType::NEGY:
			[[fallthrough]];
		case DimensionUtils::DimensionType::NEGZ:
			returnRegionEnd[dimensionLS1] =
				std::min(returnRegionBegin[dimensionLS1] + regionWidth, givenRegionEnd[dimensionLS1]);
			break;

		default:
			Log::global_log->error() << "DimensionType::ERROR received in RegionUtils::getInnerRegionSlab" << std::endl;
			mardyn_exit(1);
	}
	return {returnRegionBegin, returnRegionEnd};
}

bool RegionUtils::isMoleculeLeaving(const Molecule &molecule, const std::array<double, 3> &regionBegin,
									const std::array<double, 3> &regionEnd, DimensionUtils::DimensionType dimension,
									double timestepLength, double nextStepVelAdjustment) {
	const int ls1dim = convertEnumToLS1DimIndex(dimension);
	const int direction = findSign(dimension);
	const double newPos = molecule.r(ls1dim) + (timestepLength * (molecule.v(ls1dim) + nextStepVelAdjustment));
	if (newPos <= regionBegin[ls1dim] && direction < 0)
		return true;
	if (newPos >= regionEnd[ls1dim] && direction > 0)
		return true;
	return false;
}

std::tuple<std::array<double, 3>, std::array<double, 3>> RegionUtils::getOuterRegionSlab(
	const std::array<double, 3> &givenRegionBegin, const std::array<double, 3> &givenRegionEnd,
	DimensionUtils::DimensionType dimension, const std::array<double, 3> &regionWidth) {
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

	switch (dimension) {
		// in positive case, move the box begin to edge of domain, and box end to
		// beyond
		case DimensionUtils::DimensionType::POSX:
			[[fallthrough]];
		case DimensionUtils::DimensionType::POSY:
			[[fallthrough]];
		case DimensionUtils::DimensionType::POSZ:
			returnRegionBegin[dimensionLS1] = returnRegionEnd[dimensionLS1];
			returnRegionEnd[dimensionLS1] += regionWidth[dimensionLS1];
			break;

		// in negative case, move the box end to edge of domain, and box begin to
		// beyond
		case DimensionUtils::DimensionType::NEGX:
			[[fallthrough]];
		case DimensionUtils::DimensionType::NEGY:
			[[fallthrough]];
		case DimensionUtils::DimensionType::NEGZ:
			returnRegionEnd[dimensionLS1] = returnRegionBegin[dimensionLS1];
			returnRegionBegin[dimensionLS1] -= regionWidth[dimensionLS1];
			break;

		default:
			Log::global_log->error() << "DimensionType::ERROR received in RegionUtils::getOuterRegionSlab" << std::endl;
			mardyn_exit(1);
	}
	return std::make_tuple(returnRegionBegin, returnRegionEnd);
}
