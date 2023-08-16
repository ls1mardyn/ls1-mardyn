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

BoundaryHandler::BoundaryHandler() : _boundaries{
	{DimensionType::POSX, BoundaryType::PERIODIC}, {DimensionType::POSY, BoundaryType::PERIODIC}, 
	{DimensionType::POSZ, BoundaryType::PERIODIC}, {DimensionType::NEGX, BoundaryType::PERIODIC}, 
	{DimensionType::NEGY, BoundaryType::PERIODIC}, {DimensionType::NEGZ, BoundaryType::PERIODIC},
	{DimensionType::ERROR, BoundaryType::ERROR}
	} {}


BoundaryType BoundaryHandler::getGlobalWall(DimensionType dimension) const
{
	return _boundaries.at(dimension);
}

BoundaryType BoundaryHandler::getGlobalWall(std::string dimension) const
{
	DimensionType convertedDimension = BoundaryUtils::convertStringToDimension(dimension);
	return getGlobalWall(convertedDimension);
}

BoundaryType BoundaryHandler::getGlobalWall(int dimension) const 
{
	return getGlobalWall(BoundaryUtils::convertLS1DimsToDimensionPos(dimension));
}

void BoundaryHandler::setGlobalWall(DimensionType dimension, BoundaryType value)
{
	if(dimension != DimensionType::ERROR)
		_boundaries[dimension] = value;
}

void BoundaryHandler::setGlobalWall(std::string dimension, BoundaryType value) 
{
	DimensionType convertedDimension = BoundaryUtils::convertStringToDimension(dimension);
	setGlobalWall(convertedDimension, value);
}

void BoundaryHandler::setGlobalRegion(double* start, double* end) 
{
	for(short int i = 0; i < 3; i++) {
		_globalRegionStart[i] = start[i];
		_globalRegionEnd[i] = end[i];
	}
}

void BoundaryHandler::setLocalRegion(double* start, double* end) 
{
	for(short int i = 0; i < 3; i++) {
		_localRegionStart[i] = start[i];
		_localRegionEnd[i] = end[i];
	}
}

void BoundaryHandler::setGlobalRegion(std::array<double, 3> start, std::array<double, 3> end) {
	_globalRegionStart = start;
	_globalRegionEnd = end;
}

void BoundaryHandler::setLocalRegion(std::array<double, 3> start, std::array<double, 3> end) {
	_localRegionStart = start;
	_localRegionEnd = end;
}

void BoundaryHandler::findOuterWallsInLocalRegion() 
{
	_isOuterWall[DimensionType::POSX] = (_localRegionEnd[0] == _globalRegionEnd[0]);
	_isOuterWall[DimensionType::NEGX] = (_localRegionStart[0] == _globalRegionStart[0]);
	_isOuterWall[DimensionType::POSY] = (_localRegionEnd[1] == _globalRegionEnd[1]);
	_isOuterWall[DimensionType::NEGY] = (_localRegionStart[1] == _globalRegionStart[1]);
	_isOuterWall[DimensionType::POSZ] = (_localRegionEnd[2] == _globalRegionEnd[2]);
	_isOuterWall[DimensionType::NEGZ] = (_localRegionStart[2] == _globalRegionStart[2]);
}

bool BoundaryHandler::hasInvalidBoundary() const
{
	return _boundaries.at(DimensionType::POSX) == BoundaryType::ERROR ||
		_boundaries.at(DimensionType::POSY) == BoundaryType::ERROR ||
		_boundaries.at(DimensionType::POSZ) == BoundaryType::ERROR ||
		_boundaries.at(DimensionType::NEGX) == BoundaryType::ERROR ||
		_boundaries.at(DimensionType::NEGY) == BoundaryType::ERROR ||
		_boundaries.at(DimensionType::NEGZ) == BoundaryType::ERROR;
}

bool BoundaryHandler::isOuterWall(DimensionType dimension) const
{
	return _isOuterWall.at(dimension);
}

bool BoundaryHandler::isOuterWall(int dimension) const
{
	return isOuterWall(BoundaryUtils::convertLS1DimsToDimensionPos(dimension));
}

bool BoundaryHandler::processOuterWallLeavingParticles()
{
	auto moleculeContainer = global_simulation->getMoleculeContainer(); // :-(
	double timestepLength = (global_simulation->getIntegrator())->getTimestepLength(); 
	double cutoff = moleculeContainer->getCutoff();
	for (auto const& currentWall : _isOuterWall)
	{
		if(!currentWall.second) // not an outer wall
			continue;

		switch(getGlobalWall(currentWall.first))
		{
			case BoundaryType::PERIODIC:
				//default behaviour
				break;
			
			case BoundaryType::OUTFLOW:
			case BoundaryType::REFLECTING:
			{
				//create region
				std::array<double,3> curWallRegionBegin, curWallRegionEnd;
				std::tie(curWallRegionBegin, curWallRegionEnd) = BoundaryUtils::getInnerBuffer(_localRegionStart, _localRegionEnd, currentWall.first, cutoff);
				//conversion
				const double cstylerbegin[] = {curWallRegionBegin[0], curWallRegionBegin[1], curWallRegionBegin[2]}; 
				const double cstylerend[] = {curWallRegionEnd[0], curWallRegionEnd[1], curWallRegionEnd[2]};
				
				auto particlesInRegion = moleculeContainer->regionIterator(cstylerbegin, cstylerend, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				for (auto it = particlesInRegion; it.isValid(); ++it)
				{
					Molecule curMolecule = *it;
					
					int currentDim = BoundaryUtils::convertDimensionToLS1Dims(currentWall.first);
					double halfTimestep = .5 * timestepLength;
					double halfTimestepByMass = halfTimestep / it->mass();
					double force = it->F(currentDim);
					double nextStepVelAdjustment = halfTimestepByMass * force;

					if (BoundaryUtils::isMoleculeLeaving(curMolecule, curWallRegionBegin, curWallRegionEnd, currentWall.first, timestepLength, nextStepVelAdjustment))
					{
						if(getGlobalWall(currentWall.first) == BoundaryType::REFLECTING)
						{
							double currentVel = it->v(currentDim);
							it->setv(currentDim, -currentVel-nextStepVelAdjustment-nextStepVelAdjustment);
						}
						else // outflow
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

void BoundaryHandler::removeNonPeriodicHalos()
{
	auto moleculeContainer = global_simulation->getMoleculeContainer();
	double buffers[] = {moleculeContainer->get_halo_L(0) + moleculeContainer->getSkin(), 
						moleculeContainer->get_halo_L(1) + moleculeContainer->getSkin(), 
						moleculeContainer->get_halo_L(2) + moleculeContainer->getSkin()};
	for (auto const& currentWall : _isOuterWall)
	{
		if(!currentWall.second) //not an outer wall
			continue;

		switch(getGlobalWall(currentWall.first))
		{
			case BoundaryType::PERIODIC:
				//default behaviour
				break;
			
			case BoundaryType::OUTFLOW:
			case BoundaryType::REFLECTING:
			{
				//create region
				std::array<double,3> curWallRegionBegin, curWallRegionEnd;
				std::tie(curWallRegionBegin, curWallRegionEnd) = BoundaryUtils::getOuterBuffer(_localRegionStart, _localRegionEnd, currentWall.first, buffers);
				//conversion
				const double cstylerbegin[] = {curWallRegionBegin[0], curWallRegionBegin[1], curWallRegionBegin[2]}; 
				const double cstylerend[] = {curWallRegionEnd[0], curWallRegionEnd[1], curWallRegionEnd[2]};

				auto particlesInRegion = moleculeContainer->regionIterator(cstylerbegin, cstylerend, ParticleIterator::ALL_CELLS);
				for (auto it = particlesInRegion; it.isValid(); ++it)
				{
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