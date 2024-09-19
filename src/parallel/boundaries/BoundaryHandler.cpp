/*
 * BoundaryHandler.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include "BoundaryHandler.h"

#include "integrators/Integrator.h"

#include "utils/Logger.h"
#include "utils/Math.h"
#include "utils/mardyn_assert.h"

#include <algorithm>

BoundaryUtils::BoundaryType BoundaryHandler::getGlobalWallType(DimensionUtils::DimensionType dimension) const {
	return _boundaries.at(dimension);
}

BoundaryUtils::BoundaryType BoundaryHandler::getGlobalWallType(int dimension) const {
	return getGlobalWallType(DimensionUtils::convertLS1DimIndexToEnumPositive(dimension));
}

void BoundaryHandler::setGlobalWallType(DimensionUtils::DimensionType dimension, BoundaryUtils::BoundaryType value) {
	if (dimension != DimensionUtils::DimensionType::ERROR)
		_boundaries[dimension] = value;
	else {
		Log::global_log->error() << "DimensionType::ERROR received in BoundaryHandler::setGlobalWallType!!"
								 << std::endl;
		mardyn_exit(1);
	}
}

void BoundaryHandler::setGlobalRegion(const double *start, const double *end) {
	for (short int i = 0; i < 3; i++) {
		_globalRegionStart[i] = start[i];
		_globalRegionEnd[i] = end[i];
	}
}

void BoundaryHandler::setLocalRegion(const double *start, const double *end) {
	for (short int i = 0; i < 3; i++) {
		_localRegionStart[i] = start[i];
		_localRegionEnd[i] = end[i];
	}
}

void BoundaryHandler::updateGlobalWallLookupTable() {
	_isGlobalWall[DimensionUtils::DimensionType::POSX] = isNearRel(_localRegionEnd[0], _globalRegionEnd[0]);
	_isGlobalWall[DimensionUtils::DimensionType::NEGX] = isNearRel(_localRegionStart[0], _globalRegionStart[0]);
	_isGlobalWall[DimensionUtils::DimensionType::POSY] = isNearRel(_localRegionEnd[1], _globalRegionEnd[1]);
	_isGlobalWall[DimensionUtils::DimensionType::NEGY] = isNearRel(_localRegionStart[1], _globalRegionStart[1]);
	_isGlobalWall[DimensionUtils::DimensionType::POSZ] = isNearRel(_localRegionEnd[2], _globalRegionEnd[2]);
	_isGlobalWall[DimensionUtils::DimensionType::NEGZ] = isNearRel(_localRegionStart[2], _globalRegionStart[2]);
}

bool BoundaryHandler::hasGlobalInvalidBoundary() const {
	return std::any_of(_boundaries.begin(), _boundaries.end(), [](const auto &keyVal) {
		const auto [dim, boundaryType] = keyVal;
		return boundaryType == BoundaryUtils::BoundaryType::ERROR;
	});
}

bool BoundaryHandler::hasGlobalNonPeriodicBoundary() const {
	return std::any_of(_boundaries.begin(), _boundaries.end(), [](const auto &keyVal) {
		const auto [dim, boundaryType] = keyVal;
		return boundaryType == BoundaryUtils::BoundaryType::PERIODIC;
	});
}

bool BoundaryHandler::isGlobalWall(DimensionUtils::DimensionType dimension) const {
	return _isGlobalWall.at(dimension);
}

bool BoundaryHandler::isGlobalWall(int dimension) const {
	return isGlobalWall(DimensionUtils::convertLS1DimIndexToEnumPositive(dimension));
}

void BoundaryHandler::processGlobalWallLeavingParticles(ParticleContainer *moleculeContainer,
														double timestepLength) const {
	const auto cutoff = moleculeContainer->getCutoff();
	for (auto const [currentDim, currentWallIsGlobalWall] : _isGlobalWall) {
		if (!currentWallIsGlobalWall)
			continue;

		switch (getGlobalWallType(currentDim)) {
		case BoundaryUtils::BoundaryType::PERIODIC:
			// nothing changes from normal ls1 behaviour, so leaving particles not touched by BoundaryHandler and are
			// processed by DomainDecompBase::handleDomainLeavingParticles()
			break;

		case BoundaryUtils::BoundaryType::OUTFLOW:
			[[fallthrough]];
		case BoundaryUtils::BoundaryType::REFLECTING: {
			// create region by using getInnerRegionSlab()
			const auto [curWallRegionBegin, curWallRegionEnd] =
				RegionUtils::getInnerRegionSlab(_localRegionStart, _localRegionEnd, currentDim, cutoff);
			// grab an iterator from the converted coords
			const auto particlesInRegion = moleculeContainer->regionIterator(
				curWallRegionBegin.data(), curWallRegionEnd.data(), ParticleIterator::ONLY_INNER_AND_BOUNDARY);

			// iterate through all molecules
			for (auto moleculeIter = particlesInRegion; moleculeIter.isValid(); ++moleculeIter) {
				// Calculate the change in velocity, which the leapfrog method will
				// apply in the next velocity update to the dimension of interest.
				const int currentDimInt = DimensionUtils::convertEnumToLS1DimIndex(currentDim);
				const double halfTimestep = .5 * timestepLength;
				const double halfTimestepByMass = halfTimestep / moleculeIter->mass();
				const double force = moleculeIter->F(currentDimInt);
				const double nextStepVelAdjustment = halfTimestepByMass * force;

				// check if the molecule would leave the bounds
				if (RegionUtils::isMoleculeLeaving(*moleculeIter, curWallRegionBegin, curWallRegionEnd, currentDim,
												   timestepLength, nextStepVelAdjustment)) {
					if (getGlobalWallType(currentDim) == BoundaryUtils::BoundaryType::REFLECTING) {
						const double currentVel = moleculeIter->v(currentDimInt);
						// change the velocity in the dimension of interest such that when
						// the leapfrog integrator adds nextStepVelAdjustment in the next
						// velocity update, the final result ends up being the intended,
						// reversed velocity: -(currentVel+nextStepVelAdjustment)
						moleculeIter->setv(currentDimInt, -currentVel - nextStepVelAdjustment - nextStepVelAdjustment);
					} else { // outflow, delete the particle if it would leave
						moleculeContainer->deleteMolecule(moleculeIter, false);
					}
				}
			}
			break;
		}
		default:
			Log::global_log->error()
				<< "BoundaryType::ERROR received in BoundaryHandler::processGlobalWallLeavingParticles!" << std::endl;
			mardyn_exit(1);
		}
	}
}

void BoundaryHandler::removeNonPeriodicHalos(ParticleContainer *moleculeContainer) const {
	// get halo lengths in each dimension
	const std::array<double, 3> haloWidths = {moleculeContainer->getHaloWidthForDimension(0),
											  moleculeContainer->getHaloWidthForDimension(1),
											  moleculeContainer->getHaloWidthForDimension(2)};
	for (auto const [currentDim, currentWallIsGlobalWall] : _isGlobalWall) {
		if (!currentWallIsGlobalWall)
			continue;

		switch (getGlobalWallType(currentDim)) {
		case BoundaryUtils::BoundaryType::PERIODIC:
			// nothing changes from normal ls1 behaviour, so empty case, and halo particles left untouched
			break;

		case BoundaryUtils::BoundaryType::OUTFLOW:
			[[fallthrough]];
		case BoundaryUtils::BoundaryType::REFLECTING: {
			// create region by using getOuterRegionSlab()
			auto const [curWallRegionBegin, curWallRegionEnd] =
				RegionUtils::getOuterRegionSlab(_localRegionStart, _localRegionEnd, currentDim, haloWidths);

			// grab an iterator from the converted coords
			auto particlesInRegion = moleculeContainer->regionIterator(
				curWallRegionBegin.data(), curWallRegionEnd.data(), ParticleIterator::ALL_CELLS);
			for (auto moleculeIter = particlesInRegion; moleculeIter.isValid(); ++moleculeIter) {
				// delete all halo particles
				moleculeContainer->deleteMolecule(moleculeIter, false);
			}
			break;
		}
		default:
			Log::global_log->error() << "BoundaryType::ERROR received in BoundaryHandler::removeNonPeriodicHalos!"
									 << std::endl;
			mardyn_exit(1);
		}
	}
}
