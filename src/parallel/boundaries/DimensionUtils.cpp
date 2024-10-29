/*
 * DimensionUtils.cpp
 *
 *  Created on: 18 Sep 2024
 *      Author: amartyads
 */

#include "DimensionUtils.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h" // for MARDYN_EXIT()

#include <sstream> // for ostringstream

bool DimensionUtils::isDimensionNumericPermissible(int dim) {
	// permissible dimensions are {-1, -2, -3, 1, 2, 3}
	return (dim >= -3 && dim <= 3 && dim != 0);
}

DimensionUtils::DimensionType DimensionUtils::convertNumericToDimension(int dim) {
	if (!isDimensionNumericPermissible(dim)) {
		std::ostringstream error_message;
		error_message << "Invalid dimension passed for enum conversion. Received value: " << dim
								 << std::endl;
		MARDYN_EXIT(error_message.str());
		return DimensionType::ERROR;
	}
	switch (dim) {
		case 1:
			return DimensionType::POSX;
		case 2:
			return DimensionType::POSY;
		case 3:
			return DimensionType::POSZ;
		case -1:
			return DimensionType::NEGX;
		case -2:
			return DimensionType::NEGY;
		case -3:
			return DimensionType::NEGZ;
	}
	return DimensionType::ERROR; // warning suppression
}

DimensionUtils::DimensionType DimensionUtils::convertLS1DimIndexToEnumPositive(int dim) {
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

std::string DimensionUtils::convertDimensionToString(DimensionType dimension) {
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
		default: // ERROR
			std::ostringstream error_message;
			error_message << "DimesionType::ERROR received in DimensionUtils::convertDimensionToString!"
									 << std::endl;
			MARDYN_EXIT(error_message.str());
			return "error";
	}
}

std::string DimensionUtils::convertDimensionToStringAbs(DimensionType dimension) {
	return convertDimensionToString(dimension).substr(1, 1);
}

int DimensionUtils::convertDimensionToNumeric(DimensionType dimension) {
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
			std::ostringstream error_message;
			error_message << "DimesionType::ERROR received in DimensionUtils::convertDimensionToNumeric!"
									 << std::endl;
			MARDYN_EXIT(error_message.str());
			return 0;
	}
}

int DimensionUtils::convertDimensionToNumericAbs(DimensionType dimension) {
	return std::abs(convertDimensionToNumeric(dimension));
}

int DimensionUtils::convertEnumToLS1DimIndex(DimensionType dimension) {
	return convertDimensionToNumericAbs(dimension) - 1;
}
