/**
 * @file GeneralDomainDecomposition.h
 * @author seckler, Georg von Bismarck
 * @date 27.10.2026
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

	// documentation see base class (DomainDecompBase.h)
	~GeneralDomainDecomposition() override;

	// Read in XML configuration for GeneralDomainDecomposition and all its included objects.
	void readXML(XMLfileUnits& xmlconfig) override;

	// documentation see base class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain) override;

	// documentation see base class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain) override;

	void balanceAndExchange(double lastTraversalTime, bool forceRebalancing, ParticleContainer* moleculeContainer,
							Domain* domain) override;

	// returns a vector of the neighbour ranks in x y and z direction (only neighbours connected by an area to local
	// area)
	std::vector<int> getNeighbourRanks() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanks() not yet implemented");
	}

	// documentation see base class (DomainDecompBase.h)
	std::vector<int> getNeighbourRanksFullShell() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanksFullShell() not yet implemented");
	}

	// documentation see base class (DomainDecompBase.h)
	void prepareNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								 unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation see base class (DomainDecompBase.h)
	void finishNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation see base class (DomainDecompBase.h)
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
	 * Return the coefficients of variation (CV) also known as relative standard deviation (RSD)
	 * The coefficient of variation is calculated by dividing the standard deviation by the mean
	 * @param data array of double
	 * @param size size of array
	 * @return coefficients of variation 
	 */
	double getCV(double* data, const int size);
	
	/** 
	 * Return the Max divided Min of an array
	 * @param data array of double
	 * @param size size of array
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
	 * Check, based on the specified settings, whether rebalancing should be carried out in the specified iteration
	 * @param step current step of the simulation
	 */
	bool checkRebalancing(size_t step);
	
	/**
	 * Checks, based on available parameters (e.g. lastTraversalTime), whether rebalancing appears to be sensible. 
	 * This function uses the _imbalanceThresholdMode value internally to decide which method to use for the assessment. 
	 * This value, along with other values required for the assessment, is set via the settings.  
	 * If the _imbalanceThresholdMode value is 0, True is always returned. 
	 * For value 1, the coefficient of variation of all lastTraversalTime values is compared with the threshold (_imbalanceThresholdCV).
	 *  -> Using the getCV helper function
	 *  -> Note: The coefficient of variation was preferred over the standard deviation, as it is relative to the lastTraversalTime and is therefore a unitless quantity. 
	 * For value 2, the maximum divided by the minimum of all lastTraversalTime values is compared with the threshold (_imbalanceThresholdMinMax).
	 *  -> Using the getMaxdivMin helper function
	 * @param lastTraversalTime
	 * @return True if rebalancing appears sensible or if no method of classification is enabled. 
	 */
	bool checkNeedRebalance(double lastTraversalTime);

	/**
	 * Checks whether a proposed new domain decomposition is sensible. 
	 * To this end, the function compares the proposed decomposition with the current and previous ones. 
	 *  -> Using the domainDecompositionPercentageOfRepeatedChanges helper function
	 * The function builds a global view of the domain decomposition for the specified objective.
	 * @param proposedLocalDomain new proposed subdomain
	 * @return True if Either Rebalance is sensible or the Featcher is disabled.
	 */
	bool checkForSensibleRebalance(const DomainBox& proposedLocalDomain);

	/**
	 * Can modify the proposed new domain decomposition according to user-defined methods.
	 * Currently, only the functionality provided by the latchToGridSize function is supported.
	 * @param proposedLocalDomain new proposed subdomain
	 * @return possibly modified proposed subdomain
	 */
	DomainBox reviseNewRebalance(DomainBox proposedLocalDomain);

	/**
	 * checked whether the data in _minimalDomainSize is valid
	 * @param minimalDomainBoundary minimal DomainSize in each dimension
	 */
	void checkMinimalDomainSize(double minimalDomainBoundary);
	
	/**
	 * Calculates the geometric similarity between the previous domain decomposition change and the current proposed domain decomposition change.
	 * Geometric similarity is a procedural similarity measure defined as the geometric intersection divided by the geometric union. 
	 * In this specific context, a value of 1 would mean that the proposed change is an undo of the preceding change. 
	 *
	 * The function first generates a geometric representation of the domain decomposition change by finding, for two consecutive domain decompositions, all intersections between two subdomains from different processes. 
	 * The geometries thus identified correspond to regions that were swapped during load balancing. 
	 * This is done, firstly, between the previous and current domain decompositions, and secondly, between the current and proposed domain decompositions.
	 * The geometric similarity is then calculated between these two domain decomposition changes
	 */
	double domainDecompositionPercentageOfRepeatedChanges();

    // ------ private data members ------
	bool _cartCommunicatorCreated = false; // Indicates whether a communicator with topology information has already been created.
	
	int _imbalanceThresholdMode{0}; // 0 == disabled, 1 == Coefficient of variation (CV), 2 == MinMax 
	double _imbalanceThresholdCV{0};
	double _imbalanceThresholdMinMax{0};
	
	// Variables for storing the global view of the domain decompositions. 
	std::vector<DomainBox> _previousDomainDecomposition;
    std::vector<DomainBox> _currentDomainDecomposition;
	std::vector<DomainBox> _proposedDomainDecomposition;

	// These vectors are retained as private members for the
	// maximumRepeatedLoadChange feature to prevent repeated allocations.
	std::vector<double> _mpiBufferDomainDecomposition;
	std::vector<DomainBox> _previousDomainDecompositionChange; // between previous and current
    std::vector<DomainBox> _proposedDomainDecompositionChange; // between future and current
    
	double _maximumRepeatedLoadChange{1}; // represents a percentage, 1 == disabled

	DomainBox _localDomain{}; // The currently used local domain
	DomainBox _proposedLocalDomain{}; // The proposed version of the local domain used, without modifications

	//skin is exending minimum domain size for _minimalDomainSize
	double _skin; //is passed upon instantiation; defined by ParticleContainer settings
	DomainPoint _skinDims;  
	std::vector<double> _minimalDomainSize = {0., 0., 0.};
	
	DomainPoint _domainLength;
	double _cutoffRadius;
	
	// Setting for the funtion checkRebalancing
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
