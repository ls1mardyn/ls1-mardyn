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

void BoundaryHandler::setGlobalRegion(double* start, double* end) 
{
	for(short int i = 0; i < 3; i++) {
		globalRegionStart[i] = start[i];
		globalRegionEnd[i] = end[i];
	}
}

void BoundaryHandler::setLocalRegion(double* start, double* end) 
{
	for(short int i = 0; i < 3; i++) {
		localRegionStart[i] = start[i];
		localRegionEnd[i] = end[i];
	}
}

void BoundaryHandler::setGlobalRegion(std::array<double, 3> start, std::array<double, 3> end) {
	globalRegionStart = start;
	globalRegionEnd = end;
}

void BoundaryHandler::setLocalRegion(std::array<double, 3> start, std::array<double, 3> end) {
	localRegionStart = start;
	localRegionEnd = end;
}

void BoundaryHandler::findBoundariesInLocalRegion() 
{
	isOuterWall[DimensionType::POSX] = (localRegionEnd[0] == globalRegionEnd[0]);
	isOuterWall[DimensionType::NEGX] = (localRegionStart[0] == globalRegionStart[0]);
	isOuterWall[DimensionType::POSY] = (localRegionEnd[1] == globalRegionEnd[1]);
	isOuterWall[DimensionType::NEGY] = (localRegionStart[1] == globalRegionStart[1]);
	isOuterWall[DimensionType::POSZ] = (localRegionEnd[2] == globalRegionEnd[2]);
	isOuterWall[DimensionType::NEGZ] = (localRegionStart[2] == globalRegionStart[2]);
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

bool BoundaryHandler::processBoundaries()
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
				std::tie(curWallRegionBegin, curWallRegionEnd) = BoundaryUtils::getInnerBuffer(localRegionStart, localRegionEnd, currentWall.first, cutoff);
				//conversion
				const double cstylerbegin[] = {curWallRegionBegin[0], curWallRegionBegin[1], curWallRegionBegin[2]}; 
				const double cstylerend[] = {curWallRegionEnd[0], curWallRegionEnd[1], curWallRegionEnd[2]};
				
				auto particlesInRegion = moleculeContainer->regionIterator(cstylerbegin, cstylerend, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				for (auto it = particlesInRegion; it.isValid(); ++it)
				{
					Molecule curMolecule = *it;
					global_log->info() << "Boundary particle found " << std::endl;
					if (BoundaryUtils::isMoleculeLeaving(curMolecule, curWallRegionBegin, curWallRegionEnd, currentWall.first, timestepLength))
					{
						global_log->info() << "Boundary particle found leaving" << std::endl;
						int currentDim = BoundaryUtils::convertDimensionToLS1Dims(currentWall.first);
						if(getBoundary(currentWall.first) == BoundaryType::REFLECTING)
						{
							global_log->info() << "Reflection particle found " << std::endl;
							double vel = it->v(currentDim);
							it->setv(currentDim, -vel);
						}
						else
						{
							global_log->info() << "Outflow particle found " << std::endl;
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

void BoundaryHandler::removeHalos()
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
				std::tie(curWallRegionBegin, curWallRegionEnd) = BoundaryUtils::getInnerBuffer(localRegionStart, localRegionEnd, currentWall.first, cutoff);
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
	#ifndef MARDYN_AUTOPAS
	moleculeContainer->updateBoundaryAndHaloMoleculeCaches();
	#endif
}