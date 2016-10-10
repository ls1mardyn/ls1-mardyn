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

	// documentation see father class (DomainDecompBase.h)
	~DomainDecomposition();

	/** @brief Read in XML configuration for DomainDecomposition and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <parallelisation type="DomainDecomposition">
	     <!-- no parameters -->
	   </parallelisation>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	// documentation see father class (DomainDecompBase.h)
	bool procOwnsPos(double x, double y, double z, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain);

	void balanceAndExchange(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain);

	//! @brief writes information about the current decomposition into the given file
	//!
	//! This decomposition first writes a small header with the following lines:
	//! - "size", followed by three double values for the size of the domain
	//! - "cells", followed by three int values for the number of cells in each direction.
	//!            Each cell here stands for the domain of one process, so the product of
	//!            of those three int values is the number of processes
	//! - "procs", follwed by one int value (number of processes)
	//! - "data DomainDecomp": just this text to state that now the header is finished and data follows
	//! - one line per proc, containing the "cellid" using a x-y-z-ordering and the rank of the proc
	//!
	//! An example file for 8 procs could look like this:
	//!
	//! |size 62.0 62.0 62.0 \n
	//!  cells 2 2 2 \n
	//!  procs 8 \n
	//!  data DomainDecomp \n
	//!  1 1 \n
	//!  2 2 \n
	//!  3 3 \n
	//!  4 4 \n
	//!  5 5 \n
	//!  6 6 \n
	//!  7 7 \n
	//!  8 8
	//! @param filename name of the file into which the data will be written
	//! @param domain e.g. needed to get the bounding boxes
	void printDecomp(std::string filename, Domain* domain);

	void initCommunicationPartners(double cutoffRadius, Domain * domain);

    //returns a vector of the neighbour ranks in x y and z direction (only neighbours connected by an area to local area)
	std::vector<int> getNeighbourRanks();
	
    //returns a vector of all 26 neighbour ranks in x y and z direction
	std::vector<int> getNeighbourRanksFullShell();


	// documentation in base class
	virtual int getNonBlockingStageCount();

	// documentation in base class
	virtual void prepareNonBlockingStage(bool forceRebalancing,
			ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber);

	// documentation in base class
	virtual void finishNonBlockingStage(bool forceRebalancing,
			ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber);

	// documentation in base class
	virtual bool queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain);


private:

	//! Number of processes in each dimension (i.e. 2 for 8 processes)
	int _gridSize[DIMgeom];

	//! Grid coordinates of process
	int _coords[DIMgeom];
};

#endif /* DOMAINDECOMPOSITION_H_ */
