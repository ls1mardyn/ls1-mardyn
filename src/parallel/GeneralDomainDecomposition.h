/**
 * @file GeneralDomainDecomposition.h
 * @author seckler, Georg von Bismarck
 * @date 23.09.2026
 */

#pragma once

#include <optional>
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "DomainDecompMPIBase.h"
#include "LoadBalancer.h"



/**
 * This decomposition is meant to be able to call arbitrary load balancers.
 */
class GeneralDomainDecomposition : public DomainDecompMPIBase {
public:
    /**
	 * Constructor for the GeneralDomainDecomposition.
	 * @param cutoffRadius
	 * @param skin
	 * @param domain
	 */
	GeneralDomainDecomposition(double cutoffRadius, double skin, Domain* domain);

	GeneralDomainDecomposition(double cutoffRadius, double skin, Domain* domain, MPI_Comm comm);

	// documentation see father class (DomainDecompBase.h)
	~GeneralDomainDecomposition() override;

	// Read in XML configuration for GeneralDomainDecomposition and all its included objects.
	void readXML(XMLfileUnits& xmlconfig) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain) override;

	void balanceAndExchange(double lastTraversalTime, bool forceRebalancing, ParticleContainer* moleculeContainer,
							Domain* domain) override;

	// returns a vector of the neighbour ranks in x y and z direction (only neighbours connected by an area to local
	// area)
	std::vector<int> getNeighbourRanks() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanks() not yet implemented");
	}

	// documentation see father class (DomainDecompBase.h)
	std::vector<int> getNeighbourRanksFullShell() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanksFullShell() not yet implemented");
	}

	// documentation see father class (DomainDecompBase.h)
	void prepareNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								 unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation see father class (DomainDecompBase.h)
	void finishNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation see father class (DomainDecompBase.h)
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
	void initializeALLLoadBalancer();
	
	/**
	Creates a new MPI communicator with topology information added.
	*/
	void initMPIGridDims();

	/**
	initialize Domain Decomposition as regular grid!
	*/
    void initializeRegularGrid(const DomainPoint& domainLength, const DomainGridPoint& gridSize, const DomainGridPoint& gridCoords);
	
	/**
	 * Checks whether it is necessary to perform a rebalance
	 * @param lastTraversalTime
	 */
	bool checkNeedRebalance(double lastTraversalTime);
	
	/** 
	* Return the coefficients of variation 
	* @return coefficients of variation 
	*/
	double getCV(double* data, const int size);
	
	/** 
	* Return the Max divided Min
	* @return Max divided Min
	*/
	double getMaxdivMin(double* data, const int size); 

	/**
	 * Initializes communication partners
	 * @param moleculeContainer
	 * @param domain
	 */
	void initCommunicationPartners(Domain* domain, ParticleContainer* moleculeContainer);

	/**
	 * Calculate new distribution on process and migrate accordingly.
	 * @param lastTraversalTime time of last calculation
	 * @param domain
	 * @param particleContainer
	 */
	void rebalance(double lastTraversalTime, ParticleContainer* moleculeContainer, Domain* domain);

	/**
	 * Exchange the particles, s.t., particles are withing the particleContainer of the process they belong to.
	 * This function will rebuild the particleContainer.
	 * @param domain
	 * @param particleContainer
	 * @param newLocalDomain new local subdomain
	 */
	void migrateParticles(Domain* domain, ParticleContainer* particleContainer, DomainBox newLocalDomain);

	/**
	 * Check whether a rebalancing is necessary.
	 * @param step current step of the simulation
	 */
	bool checkRebalancing(size_t step);
	
	/**
	 * checked whether the reallocation is sensible according to specified characteristics
	 * @param proposedLocalDomain new proposed subdomain
	 * @return True if Either Rebalance is sensible or the Featcher is disabled.
	 */
	bool checkForSensibleRebalance(const DomainBox& proposedLocalDomain);

	/**
	 * TODO 
	 */
	DomainBox reviseNewRebalance(DomainBox proposedLocalDomain);

	/**
	 * checked whether the data in _minimalDomainSize is valid
	 * @param minimalDomainBoundary minimal DomainSize in each dimension
	 */
	void checkMinimalDomainSize(double minimalDomainBoundary);
	
	/**
	 * Latches domain boundaries to a grid, which is defined by _latchedGridSize.
	 * If boxMax matches the top boundary, it is not changed.
	 * @param proposedLocalDomain new subdomain
	 * @return The Latched version of the subdomain.
	 */
	DomainBox latchToGridSize(DomainBox proposedLocalDomain) const;

	/**
	 * Calculate the volume of a Bbox 
	 * @param bbox
	 *
	 * @return volumen of the Bbox
	 */
    inline double bboxVolume(const DomainBox& bbox) const;

	/**
	 * Calculates the intersection volume of two Bboxes; If there is no intersection, 0 is returned
	 * @param bbox1
	 * @param bbox2
	 *
	 * @return volumen of the intersection
	 */
	inline double bboxIntersectionVolume(const DomainBox& bbox1, const DomainBox& bbox2) const;
	
	/**
	 * Calculates the percentage of repeated changes to the total change between the previous load distribution change and the proposed one  
	*/
	double domainDecompositionPercentageOfRepeatedChanges();

    // variables
	bool _cartCommunicatorCreated = false; // Indicates whether a communicator with topology information has already been created.
	
	int _imbalanceThresholdMode{0}; // 0 == disabled, 1 == Coefficient of variation (CV), 2 == MinMax 
	double _imbalanceThresholdCV{0};
	double _imbalanceThresholdMinMax{0};
	
	// Variables for storing the total domain decompositions
	std::vector<double> _previousDomainDecomposition;
    std::vector<double> _currentDomainDecomposition;
	std::vector<double> _futureDomainDecomposition;
	// Storing the geometric differences between domain decompositions
	std::vector<DomainBox> _previousDomainDecompositionChange; // between previous and current
    std::vector<DomainBox> _futureDomainDecompositionChange; // between future and current
    
	double _maximumRepeatedLoadChange{1}; // represents a percentage, 1 == disabled

	DomainBox _localDomain{}; // The currently used local domain
	DomainBox _proposedLocalDomain{}; // The proposed version of the local domain used, without modifications

	DomainPoint _domainLength;
	std::vector<double> _minimalDomainSize = {0., 0., 0.};
	double _cutoffRadius;
	double _skin;

	size_t _steps{0};
	size_t _rebuildFrequency{10000};

	size_t _initPhase{0};
	size_t _initFrequency{500};
	
	// the LoadBalancer used
	std::unique_ptr<LoadBalancer> _loadBalancer{nullptr};
	
	// Number of processes in each dimension of the MPI process grid
	DomainGridPoint _gridSize;
	// dimensions for latching grid
	std::optional<DomainPoint> _latchedGridSize;
	
	// Coordinate of the process in the MPI process grid
	DomainGridPoint _coords;
};
