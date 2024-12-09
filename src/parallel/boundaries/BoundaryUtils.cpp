/*
 * BoundaryUtils.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include "BoundaryUtils.h"

#include "utils/Logger.h"
#include "utils/mardyn_assert.h"

#include <algorithm>
#include <cctype>
#include <sstream> // for ostringstream
#include <string>

BoundaryUtils::BoundaryType BoundaryUtils::convertStringToBoundary(const std::string &boundary) {
	std::string boundaryLowercase(boundary);
	std::transform(boundaryLowercase.begin(), boundaryLowercase.end(), boundaryLowercase.begin(),
				   [](unsigned char c) { return std::tolower(c); });
	if (boundaryLowercase.find("per") != std::string::npos)
		return BoundaryType::PERIODIC;
	if (boundaryLowercase.find("ref") != std::string::npos)
		return BoundaryType::REFLECTING;
	if (boundaryLowercase.find("out") != std::string::npos)
		return BoundaryType::OUTFLOW;
	std::ostringstream error_message;
	error_message << "Invalid boundary type passed to BoundaryUtils::convertStringToBoundary. Check your input file!"
				  << std::endl;
	MARDYN_EXIT(error_message.str());
	return BoundaryType::ERROR; // warning suppression
}

std::string BoundaryUtils::convertBoundaryToString(BoundaryType boundary) {
	switch (boundary) {
		case BoundaryType::PERIODIC:
			return "periodic";
		case BoundaryType::REFLECTING:
			return "reflecting";
		case BoundaryType::OUTFLOW:
			return "outflow";
		default:
			std::ostringstream error_message;
			error_message << "BoundaryType::ERROR received in BoundaryUtils::convertBoundaryToString!" << std::endl;
			MARDYN_EXIT(error_message.str());
	}
	return "error"; // warning suppression
}
