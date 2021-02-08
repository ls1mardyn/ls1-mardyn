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
	 * Constructor for the GeneralDomainDecomposition.
	 * @param interactionLength
	 * @param domain
	 */
	GeneralDomainDecomposition(double interactionLength, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	~GeneralDomainDecomposition() override;

	/** @brief Read in XML configuration for DomainDecomposition and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <parallelisation type="GeneralDomainDecomposition">
		  <updateFrequency>INTEGER</updateFrequency>
		  <initialPhaseTime>INTEGER</initialPhaseTime><!--time for initial rebalancing phase-->
		  <initialPhaseFrequency>INTEGER</initialPhaseFrequency><!--frequency for initial rebalancing phase-->
		  <gridSize>STRING</gridSize><!--default: 0; if non-zero, the process boundaries are fixed to multiples of
				gridSize. Comma separated string to define three different grid sizes for the different dimensions is
				possible.-->
		  <loadBalancer type="STRING"><!--STRING...type of the load balancer, currently supported: ALL-->
			<!--options for the load balancer-->
			<!--for detailed information see the readXML functions from ALLLoadBalancer.-->
		  </loadBalancer>
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
	void printDecomp(const std::string& filename, Domain* domain) override;

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
	 * Method that initializes the ALLLoadBalancer
	 */
	void initializeALL();

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
	 * Check whether a rebalancing is necessary.
	 * @param step the step number
	 * @param updateFrequency the frequency defining how often
	 * @param lastTraversalTime the time of the last traversal for this node
	 * @param initPhase Number of time steps in the initial rebalancing phase.
	 * @param initUpdateFrequency Update frequency within the initial rebalancing phase, normally higher frequency,
	 * i.e., lower value than updateFrequency.
	 * @return true if a rebuild is necessary
	 */
	static bool queryRebalancing(size_t step, size_t updateFrequency, size_t initPhase, size_t initUpdateFrequency,
								 double lastTraversalTime);

	/**
	 * Initializes communication partners
	 * @param moleculeContainer
	 * @param domain
	 */
	void initCommPartners(ParticleContainer* moleculeContainer, Domain* domain);

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

	/**
	 * Latches domain boundaries (given as boxMin and boxMax) to a grid, which is defined by _gridSize.
	 * If boxMax matches the top boundary, it is not changed.
	 * @param boxMin
	 * @param boxMax
	 * @return The new boundaries.
	 */
	std::pair<std::array<double, 3>, std::array<double, 3>> latchToGridSize(std::array<double, 3> boxMin,
																			std::array<double, 3> boxMax) {
		for (size_t ind = 0; ind < 3; ++ind) {
			double currentGridSize = (*_gridSize)[ind];
			// For boxmin, the lower domain boundary is 0, so that's always fine!
			boxMin[ind] = std::round(boxMin[ind] / currentGridSize) * currentGridSize;
			// update boxmax only if it isn't at the very top of the domain!
			if (boxMax[ind] != _domainLength[ind]) {
				boxMax[ind] = std::round(boxMax[ind] / currentGridSize) * currentGridSize;
			}
		}
		return {boxMin, boxMax};
	}

	// variables
	std::array<double, 3> _boxMin;
	std::array<double, 3> _boxMax;

	std::array<double, 3> _domainLength;
	double _interactionLength;

	size_t _steps{0};
	size_t _rebuildFrequency{10000};

	size_t _initPhase{0};
	size_t _initFrequency{500};

	/**
	 * Optionally safe a given grid size on which the process boundaries are bound/latched.
	 * If no value is given, it is not used.
	 */
	std::optional<std::array<double, 3>> _gridSize{};

	std::unique_ptr<LoadBalancer> _loadBalancer{nullptr};

	friend class GeneralDomainDecompositionTest;

};
