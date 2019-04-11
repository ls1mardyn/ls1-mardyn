/**
 * @file MultiSectionMethod.h
 * @author seckler
 * @date 11.04.19
 */

#pragma once

#include "DomainDecompMPIBase.h"

class MultiSectionMethod : public DomainDecompMPIBase {
public:
	/**
	 * Constructor for the MultiSectionMethod
	 * @param cutoffRadius
	 * @param domain
	 */
	MultiSectionMethod(double cutoffRadius, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	~MultiSectionMethod() override;

	/** @brief Read in XML configuration for DomainDecomposition and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <parallelisation type="MultiSectionMethod">

	   </parallelisation>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain) override;

	void balanceAndExchange(double lastTraversalTime, bool forceRebalancing, ParticleContainer* moleculeContainer,
							Domain* domain) override;

	//! @param filename name of the file into which the data will be written
	//! @param domain e.g. needed to get the bounding boxes
	void printDecomp(std::string filename, Domain* domain) override {
		throw std::runtime_error("MultiSectionMethod::printDecomp() not yet implemented");
	}

	// returns a vector of the neighbour ranks in x y and z direction (only neighbours connected by an area to local
	// area)
	std::vector<int> getNeighbourRanks() override {
		throw std::runtime_error("MultiSectionMethod::getNeighbourRanks() not yet implemented");
	}

	// returns a vector of all 26 neighbour ranks in x y and z direction
	std::vector<int> getNeighbourRanksFullShell() override {
		throw std::runtime_error("MultiSectionMethod::getNeighbourRanksFullShell() not yet implemented");
	}

	// documentation in base class
	void prepareNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								 unsigned int stageNumber) override {
		throw std::runtime_error("MultiSectionMethod::prepareNonBlockingStage() not yet implemented");
	}

	// documentation in base class
	void finishNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								unsigned int stageNumber) override {
		throw std::runtime_error("MultiSectionMethod::prepareNonBlockingStage() not yet implemented");
	}

	// documentation in base class
	bool queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
											double etime) override {
		throw std::runtime_error("MultiSectionMethod::queryBalanceAndExchangeNonBlocking() not yet implemented");
	}

	std::vector<CommunicationPartner> getNeighboursFromHaloRegion(Domain* domain, const HaloRegion& haloRegion,
																  double cutoff) override {
		throw std::runtime_error("MultiSectionMethod::getNeighboursFromHaloRegion() not yet implemented");
	}

private:

	// variables
	std::array<double, 3> _boxMin;
	std::array<double, 3> _boxMax;
	std::array<size_t, 3> _gridSize; //!< Number of processes in each dimension of the MPI process grid used by the MSM
	std::array<size_t, 3> _coords; //!< Coordinate of the process in the MPI process grid used by the MSM
	double _cutoffRadius;
};
