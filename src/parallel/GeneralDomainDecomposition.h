/**
 * @file GeneralDomainDecomposition.h
 * @author seckler
 * @date 11.04.19
 */

#pragma once

#include "DomainDecompMPIBase.h"
#include "LoadBalancer.h"

/**
 * This decomposition is meant to be able to call arbitrary load balancers.
 */
class GeneralDomainDecomposition : public DomainDecompMPIBase {
public:
	/**
	 * Constructor for the GeneralDomainDecomposition
	 * @param cutoffRadius
	 * @param domain
	 */
	GeneralDomainDecomposition(double cutoffRadius, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	~GeneralDomainDecomposition() override;

	/** @brief Read in XML configuration for DomainDecomposition and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <parallelisation type="GeneralDomainDecomposition">

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
		throw std::runtime_error("GeneralDomainDecomposition::printDecomp() not yet implemented");
	}

	// returns a vector of the neighbour ranks in x y and z direction (only neighbours connected by an area to local
	// area)
	std::vector<int> getNeighbourRanks() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanks() not yet implemented");
	}

	// returns a vector of all 26 neighbour ranks in x y and z direction
	std::vector<int> getNeighbourRanksFullShell() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanksFullShell() not yet implemented");
	}

	// documentation in base class
	void prepareNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								 unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation in base class
	void finishNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation in base class
	bool queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
											double etime) override {
		throw std::runtime_error(
			"GeneralDomainDecomposition::queryBalanceAndExchangeNonBlocking() not yet implemented");
	}

	std::vector<CommunicationPartner> getNeighboursFromHaloRegion(Domain* domain, const HaloRegion& haloRegion,
																  double cutoff) override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighboursFromHaloRegion() not yet implemented");
	}

private:
	/**
	 * Get the optimal grid for the given dimensions of the box and the number of processes.
	 * The grid is produced, s.t., the number of grid[0] * grid[1] * grid[2] == numProcs
	 * The edge lengths of the grid will resemble the lengths of the domain, i.e., the longest edge of the domain will
	 * also have the largest amount of grid points.
	 * @param domainLength
	 * @param numProcs
	 * @return
	 */
	static std::array<size_t, 3> getOptimalGrid(const std::array<double, 3>& domainLength, int numProcs);

	/**
	 * Get the coordinates from the rank.
	 * S.t. rank = x * gridSize[1]*gridSize[2] + y * gridSize[2] + z
	 * @param gridSize
	 * @param rank
	 * @return
	 */
	static std::array<size_t, 3> getCoordsFromRank(const std::array<size_t, 3>& gridSize, int rank);

	/**
	 * Returns boxMin and boxMax according to regular grid.
	 * @param domainLength
	 * @param gridSize
	 * @param gridCoords
	 * @return boxMin and boxMax
	 */
	static std::tuple<std::array<double, 3>, std::array<double, 3>> initializeRegularGrid(
		const std::array<double, 3>& domainLength, const std::array<size_t, 3>& gridSize,
		const std::array<size_t, 3>& gridCoords);

	// variables
	std::array<double, 3> _boxMin;
	std::array<double, 3> _boxMax;
	double _cutoffRadius;
	std::array<bool, 3> _coversWholeDomain{};

	size_t _steps{0};
	size_t _rebuildFrequency{10000};

	std::unique_ptr<LoadBalancer> _loadBalancer{nullptr};

	friend class GeneralDomainDecompositionTest;
};
