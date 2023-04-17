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

BoundaryType BoundaryHandler::getBoundary(int dimension) const 
{
	DimensionType toRet = DimensionType::ERROR;
	switch(dimension)
	{
		case 0:
			toRet = DimensionType::POSX;
			break;
		case 1:
			toRet = DimensionType::POSY;
			break;
		default: //case 2:
			toRet = DimensionType::POSZ;
	}
	return boundaries.at(toRet);
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

bool BoundaryHandler::processBoundaries(std::array<double,3> startRegion, std::array<double,3> endRegion)
{
	auto moleculeContainer = global_simulation->getMoleculeContainer(); // :-(
	double timestepLength = (global_simulation->getIntegrator())->getTimestepLength(); 
	double cutoff = moleculeContainer->getCutoff();
	for (auto const& currentWall : isOuterWall)
	{
		//global_log->info() << "wall number " << BoundaryUtils::convertDimensionToString(currentWall.first) << " : " << currentWall.second << std::endl;
		if(!currentWall.second)
			continue;

		switch(getBoundary(currentWall.first))
		{
			case BoundaryType::PERIODIC:
				//default behaviour
				break;
			
			case BoundaryType::OUTFLOW:
			case BoundaryType::REFLECTING:
			{
				//create region
				std::array<double,3> curWallRegionBegin, curWallRegionEnd;
				std::tie(curWallRegionBegin, curWallRegionEnd) = BoundaryUtils::getInnerBuffer(startRegion, endRegion, currentWall.first, cutoff);
				//conversion
				const double cstylerbegin[] = {curWallRegionBegin[0], curWallRegionBegin[1], curWallRegionBegin[2]}; 
				const double cstylerend[] = {curWallRegionEnd[0], curWallRegionEnd[1], curWallRegionEnd[2]};
				
				auto particlesInRegion = moleculeContainer->regionIterator(cstylerbegin, cstylerend, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				for (auto it = particlesInRegion; it.isValid(); ++it)
				{
					Molecule curMolecule = *it;
					if (BoundaryUtils::isMoleculeLeaving(curMolecule, curWallRegionBegin, curWallRegionEnd, currentWall.first, timestepLength))
					{
						int currentDim = BoundaryUtils::convertDimensionToLS1Dims(currentWall.first);
						if(getBoundary(currentWall.first) == BoundaryType::REFLECTING)
						{
							double vel = it->v(currentDim);
							it->setv(currentDim, -vel);
						}
						else
						{
							moleculeContainer->deleteMolecule(it, false);
						}
					}
				}
				break;
			}
			default:
				global_log->error() << "Boundary type error! Received type not allowed!" << std::endl;
				Simulation::exit(1);
		}
	}
	return true;
}

void BoundaryHandler::removeHalos(std::array<double,3> startRegion, std::array<double,3> endRegion)
{
	auto moleculeContainer = global_simulation->getMoleculeContainer();
	double cutoff = moleculeContainer->getCutoff();
	for (auto const& currentWall : isOuterWall)
	{
		//global_log->info() << "wall number " << BoundaryUtils::convertDimensionToString(currentWall.first) << " : " << currentWall.second << std::endl;
		if(!currentWall.second) //not an outer wall
			continue;

		switch(getBoundary(currentWall.first))
		{
			case BoundaryType::PERIODIC:
				//default behaviour
				break;
			
			case BoundaryType::OUTFLOW:
			case BoundaryType::REFLECTING:
			{
				//create region
				std::array<double,3> curWallRegionBegin, curWallRegionEnd;
				std::tie(curWallRegionBegin, curWallRegionEnd) = BoundaryUtils::getInnerBuffer(startRegion, endRegion, currentWall.first, cutoff);
				//conversion
				const double cstylerbegin[] = {curWallRegionBegin[0], curWallRegionBegin[1], curWallRegionBegin[2]}; 
				const double cstylerend[] = {curWallRegionEnd[0], curWallRegionEnd[1], curWallRegionEnd[2]};

				auto particlesInRegion = moleculeContainer->regionIterator(cstylerbegin, cstylerend, ParticleIterator::ALL_CELLS);
				for (auto it = particlesInRegion; it.isValid(); ++it)
				{
					global_log->info() << "Halo particle found " << std::endl;
					moleculeContainer->deleteMolecule(it, false);
				}
				break;
			}
			default:
				global_log->error() << "Boundary type error! Received type not allowed!" << std::endl;
				Simulation::exit(1);
		}
	}
	moleculeContainer->updateBoundaryAndHaloMoleculeCaches();
}