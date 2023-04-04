/*
 * BoundaryHandler.h
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#pragma once

#include <map>
#include <vector>
#include <string>
#include "BoundaryType.h"
#include "DimensionType.h"


class BoundaryHandler
{
public:
	BoundaryHandler();
	BoundaryType getBoundary(std::string dimension) const;
	void setBoundary(std::string dimension, BoundaryType value);
	BoundaryType getBoundary(DimensionType dimension) const;
	void setBoundary(DimensionType dimension, BoundaryType value);
	bool hasInvalidBoundary() const;

	void setOuterWalls(std::map<DimensionType, bool> isOuterWallOth) {isOuterWall = isOuterWallOth;}
	bool processBoundaries(double* startRegion, double* endRegion);

private:
	std::map<DimensionType, BoundaryType> boundaries;
	std::map<DimensionType, bool> isOuterWall;
};