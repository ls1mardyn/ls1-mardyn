#ifndef DOMAINDECOMPOSITION_H_
#define DOMAINDECOMPOSITION_H_

#include "DomainDecompMPIBase.h"


/** @brief Basic domain decomposition based parallelisation, dividing the
 * domain into \#procs equal sized cuboids
 *
 * In a domain decomposition, each process gets part of the spacial domain.
 * In this implementation, the whole domain has to be a cuboid which is
 * decomposed into several, equally sized smaller cuboids.
 * At the boundary, each process needs molecules from neighbouring domains to
 * be able to calculate the forces on the own molecules.
 * Molecules are moving across the boundaries of local domains. So methods are
 * implemented to transfer those molecules.
 *
 * @cite Griebel-2007
 */
class DomainDecomposition : public DomainDecompMPIBase {
public:
	//! @brief The constructor has to determine the own rank and the number of neighbours and
	//!        sets up the topology
	DomainDecomposition();

	DomainDecomposition(MPI_Comm comm, const std::array<int, DIMgeom> &gridSize);
	// documentation see father class (DomainDecompBase.h)
	~DomainDecomposition();

	/** @brief Read in XML configuration for DomainDecomposition and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <parallelisation type="DomainDecomposition">
	     <!-- structure handled by DomainDecompMPIBase -->
	     <MPIGridDims> <x>INT</x> <y>INT</y> <z>INT</z> </MPIGridDims>
	   </parallelisation>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain) override;

	void balanceAndExchange(double lastTraversalTime, bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain) override;

	void initCommunicationPartners(double cutoffRadius, Domain * domain,ParticleContainer* moleculeContainer);

    //returns a vector of the neighbour ranks in x y and z direction (only neighbours connected by an area to local area)
	std::vector<int> getNeighbourRanks() override;

    //returns a vector of all 26 neighbour ranks in x y and z direction
	std::vector<int> getNeighbourRanksFullShell() override;

	//returns the ranks of the neighbours
	std::vector<std::vector<std::vector<int>>> getAllRanks() override;


	// documentation in base class
	void prepareNonBlockingStage(bool forceRebalancing,
			ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber) override;

	// documentation in base class
	void finishNonBlockingStage(bool forceRebalancing,
			ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber) override;

	// documentation in base class
	bool queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain, double etime) override;

	std::vector<CommunicationPartner> getNeighboursFromHaloRegion(Domain* domain, const HaloRegion& haloRegion, double cutoff) override;

protected:
	void initMPIGridDims();

	std::array<int, DIMgeom> _gridSize; //!< Number of processes in each dimension of the MPI process grid
	int _coords[DIMgeom]; //!< Coordinate of the process in the MPI process grid
};

#endif /* DOMAINDECOMPOSITION_H_ */
