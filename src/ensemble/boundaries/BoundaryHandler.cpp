/*
 * BoundaryHandler.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include <algorithm>

#include "BoundaryHandler.h"
#include "BoundaryUtils.h"
#include "utils/Logger.h" 

BoundaryHandler::BoundaryHandler() : boundaries{
	{DimensionType::POSX, BoundaryType::PERIODIC}, {DimensionType::POSY, BoundaryType::PERIODIC}, 
	{DimensionType::POSZ, BoundaryType::PERIODIC}, {DimensionType::NEGX, BoundaryType::PERIODIC}, 
	{DimensionType::NEGY, BoundaryType::PERIODIC}, {DimensionType::NEGZ, BoundaryType::PERIODIC},
	{DimensionType::ERROR, BoundaryType::ERROR}
	} {}

BoundaryType BoundaryHandler::getBoundary(std::string dimension) const
{
	DimensionType convertedDimension = BoundaryUtils::convertStringToDimension(dimension);
	return boundaries.at(convertedDimension);
}

void BoundaryHandler::setBoundary(std::string dimension, BoundaryType value) 
{
	DimensionType convertedDimension = BoundaryUtils::convertStringToDimension(dimension);
	if(convertedDimension == DimensionType::ERROR)
		return;
	boundaries[convertedDimension] = value;
}

BoundaryType BoundaryHandler::getBoundary(DimensionType dimension) const
{
	return boundaries.at(dimension);
}

void BoundaryHandler::setBoundary(DimensionType dimension, BoundaryType value)
{
	if(dimension != DimensionType::ERROR)
		boundaries[dimension] = value;
}
