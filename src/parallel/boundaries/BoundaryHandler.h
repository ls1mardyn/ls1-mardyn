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
#include <array>
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
	BoundaryType getBoundary(int dimension) const;
	bool hasInvalidBoundary() const;

	void setOuterWalls(std::map<DimensionType, bool> isOuterWallOth) {isOuterWall = isOuterWallOth;}
	bool processBoundaries(std::array<double,3> startRegion, std::array<double,3> endRegion);
	void removeHalos(std::array<double,3> startRegion, std::array<double,3> endRegion);

private:
	std::map<DimensionType, BoundaryType> boundaries;
	std::map<DimensionType, bool> isOuterWall;
};