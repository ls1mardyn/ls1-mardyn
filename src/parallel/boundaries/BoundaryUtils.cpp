/*
 * BoundaryUtils.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include <algorithm>

#include "BoundaryUtils.h"
#include "utils/Logger.h" 

BoundaryUtils::BoundaryUtils() : boundaries{
	{DimensionType::POSX, BoundaryType::PERIODIC}, {DimensionType::POSY, BoundaryType::PERIODIC}, 
	{DimensionType::POSZ, BoundaryType::PERIODIC}, {DimensionType::NEGX, BoundaryType::PERIODIC}, 
	{DimensionType::NEGY, BoundaryType::PERIODIC}, {DimensionType::NEGZ, BoundaryType::PERIODIC}
	} {}

BoundaryType BoundaryUtils::getBoundary(std::string dimension) 
{
	DimensionType convertedDimension = convertStringToDimension(dimension);
	if (convertedDimension == DimensionType::ERROR)
		return BoundaryType::ERROR;
	return boundaries[convertedDimension];
}

void BoundaryUtils::setBoundary(std::string dimension, BoundaryType value) 
{
	DimensionType convertedDimension = convertStringToDimension(dimension);
	if (convertedDimension == DimensionType::ERROR)
		return;
	boundaries[convertedDimension] = value;
}

BoundaryType BoundaryUtils::getBoundary(DimensionType dimension) { return BoundaryType(); }

void BoundaryUtils::setBoundary(DimensionType dimension, BoundaryType value) {}

DimensionType BoundaryUtils::convertStringToDimension(std::string dimension) 
{ 
	if(std::find(permissibleDimensions.begin(),permissibleDimensions.end(),dimension) == permissibleDimensions.end())
	{
		Log::global_log->error() << "Invalid dimension passed for boundary check" << std::endl;
		return DimensionType::ERROR;
	}
	if(dimension == "+x")
		return DimensionType::POSX;
	if(dimension == "+y")
		return DimensionType::POSY;
	if(dimension == "+z")
		return DimensionType::POSZ;
	if(dimension == "-x")
		return DimensionType::NEGX;
	if(dimension == "-y")
		return DimensionType::NEGY;
	//if(dimension == "-z")
	return DimensionType::NEGZ;
}
