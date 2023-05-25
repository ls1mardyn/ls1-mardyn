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

	void setGlobalRegion(double* start, double* end);
	void setLocalRegion(double* start, double* end);

	void setGlobalRegion(std::array<double,3> start, std::array<double,3> end);
	void setLocalRegion(std::array<double,3> start, std::array<double,3> end);

	void findBoundariesInLocalRegion();
	bool processBoundaries();
	void removeHalos();

private:
	std::map<DimensionType, BoundaryType> boundaries;
	std::map<DimensionType, bool> isOuterWall;
	std::array<double,3> globalRegionStart, globalRegionEnd;
	std::array<double,3> localRegionStart, localRegionEnd;
};