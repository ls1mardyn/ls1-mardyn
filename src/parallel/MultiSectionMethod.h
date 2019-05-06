/**
 * @file MultiSectionMethod.h
 * @author seckler
 * @date 11.04.19
 */

#pragma once

#include "DomainDecompMPIBase.h"

/**
 * This method uses sth. similar to a grid of mpi-ranks and works in the following way:
 * 1. First the domain is split along the x dimension.
 * 2. Secondly in the y dimension
 * 3. Thirdly in the z dimension
 */
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

	/**
	 * Check whether a rebalancing is necessary
	 * @param step the step number
	 * @param updateFrequency the frequency defining how often
	 * @param lastTraversalTime the time of the last traversal for this node
	 * @return true if a rebuild is necessary
	 */
	static bool queryRebalancing(size_t step, size_t updateFrequency, double lastTraversalTime);

	/**
	 * Initializes communication partners
	 * @param moleculeContainer
	 * @param domain
	 */
	void initCommPartners(ParticleContainer* moleculeContainer, Domain* domain);

	/**
	 *
	 * @param lastTraversalTime time for last traversal
	 * @param particleContainer the ParticleContainer
	 * @param domain the domain
	 * @return tuple of new boxMin and boxMax
	 */
	std::tuple<std::array<double, 3>, std::array<double, 3>> doRebalancing(double lastTraversalTime,
																		   ParticleContainer* particleContainer,
																		   Domain* domain);

	/**
	 * Exchange the particles, s.t., particles are withing the particleContainer of the process they belong to.
	 * This function will rebuild the particleContainer.
	 * @param domain
	 * @param particleContainer
	 * @param newMin new minimum of the own subdomain
	 * @param newMax new maximum of the own subdomain
	 */
	void migrateParticles(Domain* domain, ParticleContainer* particleContainer, std::array<double, 3> newMin,
						  std::array<double, 3> newMax);

	// variables
	std::array<double, 3> _boxMin;
	std::array<double, 3> _boxMax;
	std::array<size_t, 3> _gridSize;  //!< Number of processes in each dimension of the MPI process grid used by the MSM
	std::array<size_t, 3> _gridCoords;  //!< Coordinate of the process in the MPI process grid used by the MSM
	double _cutoffRadius;
	size_t _step;
	size_t _updateFrequency;

	friend class MultiSectionMethodTest;
};
