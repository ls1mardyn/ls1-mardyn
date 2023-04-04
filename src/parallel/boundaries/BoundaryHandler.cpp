/*
 * BoundaryHandler.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include <algorithm>

#include "BoundaryHandler.h"
#include "BoundaryUtils.h"

#include "Simulation.h"
#include "integrators/Integrator.h"

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

bool BoundaryHandler::hasInvalidBoundary() const
{
	return boundaries.at(DimensionType::POSX) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::POSY) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::POSZ) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::NEGX) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::NEGY) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::NEGZ) == BoundaryType::ERROR;

}

bool BoundaryHandler::processBoundaries(double* startRegion, double* endRegion)
{
	auto moleculeContainer = global_simulation->getMoleculeContainer();
	double timestepLength = (global_simulation->getIntegrator())->getTimestepLength(); 
	double cutoff = moleculeContainer->getCutoff();
	for (auto const& currentWall : isOuterWall)
	{
		global_log->info() << "wall number " << BoundaryUtils::convertDimensionToString(currentWall.first) << " : " << currentWall.second << std::endl;
		if(!currentWall.second)
			continue;

		switch(getBoundary(currentWall.first))
		{
			case BoundaryType::PERIODIC:
				//default behaviour
				break;
			
			case BoundaryType::OUTFLOW:
				//delete exiting particles
				//remove from invalidparticles

				break;

			case BoundaryType::REFLECTING:
			{
				//create region
				double curWallRegionBegin[3], curWallRegionEnd[3];
				bool successRegionSet = BoundaryUtils::setRegionToParams(startRegion, endRegion, currentWall.first, cutoff, curWallRegionBegin, curWallRegionEnd);
				if (!successRegionSet)
				{
					global_log->error() << "Error while setting the region for boundary conditions " << std::endl;
					Simulation::exit(1);
				}
				auto particlesInRegion = moleculeContainer->regionIterator(curWallRegionBegin, curWallRegionEnd, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				for (auto it = particlesInRegion; it.isValid(); ++it)
				{
					Molecule curMolecule = *it;
					if (BoundaryUtils::isMoleculeLeaving(curMolecule, curWallRegionBegin, curWallRegionEnd, currentWall.first, timestepLength))
					{
						int currentDim = BoundaryUtils::convertDimensionToLS1Dims(currentWall.first);
						double vel = it->v(currentDim);
						it->setv(currentDim, -vel);
					}
				}
				break;
			}
			default:
				global_log->error() << "Boundary type error! Received type not allowed!" << std::endl;
				Simulation::exit(1);
		}
	}
}