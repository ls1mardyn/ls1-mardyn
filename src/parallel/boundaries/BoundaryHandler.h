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
	BoundaryType getGlobalWall(std::string dimension) const;
	void setGlobalWall(std::string dimension, BoundaryType value);
	BoundaryType getGlobalWall(DimensionType dimension) const;
	void setGlobalWall(DimensionType dimension, BoundaryType value);
	BoundaryType getGlobalWall(int dimension) const;
	bool hasInvalidBoundary() const;
	bool hasNonPeriodicBoundary() const;

	void setGlobalRegion(double* start, double* end);
	void setLocalRegion(double* start, double* end);

	void setGlobalRegion(std::array<double,3> start, std::array<double,3> end);
	void setLocalRegion(std::array<double,3> start, std::array<double,3> end);

	void findOuterWallsInLocalRegion();
	bool isOuterWall(DimensionType dimension) const;
	bool isOuterWall(int dimension) const;
	bool processOuterWallLeavingParticles();
	void removeNonPeriodicHalos();

private:
	std::map<DimensionType, BoundaryType> _boundaries;
	std::map<DimensionType, bool> _isOuterWall;
	std::array<double,3> _globalRegionStart, _globalRegionEnd;
	std::array<double,3> _localRegionStart, _localRegionEnd;
};