#ifndef DOMAINDECOMPOSITION_H_
#define DOMAINDECOMPOSITION_H_

#define DIM 3

#define LOWER  0
#define HIGHER 1

#include "parallel/DomainDecompBase.h"
#include "parallel/ParticleData.h"
#include "parallel/CollectiveCommunication.h"

#include <mpi.h>
#include <iostream>

//! @brief Basic parallelisation which divides the domain into #procs equal sized cuboids
//! @author Martin Buchholz
//!
//! In a domain decomposition, each process gets part of the spacial domain.
//! In this implementation, the whole domain has to be a cuboid which is
//! decomposed into several smaller cuboids.
//! The main difficulty is that at the boundary, each process needs molecules
//! from neighbouring domains to be able to calculate the forces on the own
//! molecules. Another problem is that molecules are moving across the
//! boundaries of local domains. So methods are implemented to transfer those molecules
//! @todo reference to a paper describing the %domain decomposition
class DomainDecomposition : public DomainDecompBase {
public:
	//! @brief The constructor has to determine the own rank and the number of neighbours and
	//!        sets up the topology
	DomainDecomposition();

	// documentation see father class (DomainDecompBase.h)
	~DomainDecomposition();

	//! @brief exchange molecules between processes
	//!
	//! molecules which aren't in the domain of their process any
	//! more are transferred to their neighbours. Additionally, the
	//! molecules for the halo-region are transferred. To reduce the number
	//! of neighbours a single process has to communicate with, particles
	//! that i.e. have to be moved to the lower right neighbour are
	//! moved to the right neighbour first and then from the right neighbour
	//! to the lower neighbour.
	//! @param moleculeContainer needed to get those molecules which have to be exchanged
	//! @param components when creating a new Molecule-object (from the recieved data),
	//!                   the Molecule-constructor needs this component vector
	//! @param domain is e.g. needed to get the size of the local domain
	void exchangeMolecules(ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain);

	//! @brief this decompositin does no balancing, it just exchanges the particles
	//!
	//! This domain decomposition devides the domain into equally sized smaller cuboids, therefore
	//! a balancing of the load is not possible. The method only has to ensure that the particles
	//! between the processes are exchanged, therefore exchangeMolecules is called.
	//! @param balance has no influence in this implementation
	//! @param moleculeContainer needed for calculating load and to get the particles
	//! @param components when creating a new Molecule-object (from the recieved data),
	//!                   the Molecule-constructor needs this component vector
	//! @param domain is e.g. needed to get the size of the local domain
	void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	bool procOwnsPos(double x, double y, double z, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	double guaranteedDistance(double x, double y, double z, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	unsigned long countMolecules(ParticleContainer* moleculeContainer, std::vector<unsigned long> &compCount);

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain);

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

	//! @brief append the molecule date of all processes to the file
	//!
	//! Currently, parallel IO isn't used.
	//! To ensure that not more than one process writes to the file at any time,
	//! there is a loop over all processes with a barrier in between
	//! @param filename name of the file into which the data will be written
	//! @param moleculeContainer all Particles from this container will be written to the file
	void writeMoleculesToFile(std::string filename, ParticleContainer* moleculeContainer);

	// documentation see father class (DomainDecompBase.h)
	int getRank(void) {
		return _rank;
	}

	// documentation see father class (DomainDecompBase.h)
	int getNumProcs();

	// documentation see father class (DomainDecompBase.h)
	void barrier() { MPI_CHECK( MPI_Barrier(_comm) ); }

	// documentation see father class (DomainDecompBase.h)
	double getTime();

	//! @brief returns total number of molecules
	unsigned Ndistribution(unsigned localN, float* minrnd, float* maxrnd);

	//! @brief checks identity of random number generators
	void assertIntIdentity(int IX);
	void assertDisjunctivity(TMoleculeContainer* mm);

	//##################################################################
	// The following methods with prefix "collComm" are all used
	// in the context of collective communication. Each of the methods
	// basically has to call the corresponding method from the class
	// CollectiveCommunication (or CollectiveCommDummy in the sequential
	// case). To get information about how to use this methods, read
	// the documentation of the class CollectiveCommunication and of the
	// father class of this class (DomainDecompBase.h)
	//##################################################################
	void collCommInit(int numValues) {
		_collComm.init(_comm, numValues);
	}

	void collCommFinalize() {
		_collComm.finalize();
	}

	void collCommAppendInt(int intValue) {
		_collComm.appendInt(intValue);
	}

	void collCommAppendUnsLong(unsigned long unsLongValue) {
		_collComm.appendUnsLong(unsLongValue);
	}

	void collCommAppendFloat(float floatValue) {
		_collComm.appendFloat(floatValue);
	}

	void collCommAppendDouble(double doubleValue) {
		_collComm.appendDouble(doubleValue);
	}

	void collCommAppendLongDouble(long double longDoubleValue) {
		_collComm.appendLongDouble(longDoubleValue);
	}

	int collCommGetInt() {
		return _collComm.getInt();
	}

	unsigned long collCommGetUnsLong() {
		return _collComm.getUnsLong();
	}

	float collCommGetFloat() {
		return _collComm.getInt();
	}

	double collCommGetDouble() {
		return _collComm.getDouble();
	}

	long double collCommGetLongDouble() {
		return _collComm.getLongDouble();
	}

	void collCommAllreduceSum() {
		_collComm.allreduceSum();
	}

	void collCommBroadcast() {
		_collComm.broadcast();
	}

private:
	//! determines and returns the rank of the process at the given coordinates
	int getRank(int x, int y, int z);
	//! with the given number of processes, the dimensions of the grid are calculated
	void setGridSize(int num_procs);

	//! new topology after initializing the torus
	MPI_Comm _comm;
	int _comm_size;
	MPI_Group _comm_group;
	MPI_Group _neighbours_groups[DIM][2];

	MPI_Datatype _mpi_Particle_data;
	//! Number of processes in each dimension (i.e. 2 for 8 processes)
	int _gridSize[DIM];
	//! Grid coordinates of process
	int _coords[DIM];
	//!  rank of process
	int _rank;
	//! Array of neighbour ranks.
	//! The first array index specifies the coordinate index,
	//! the second one the direction. For the later use the predefined LOWER and HIGHER macros.
	int _neighbours[DIM][2];

	//! variable used for different kinds of collective operations
	CollectiveCommunication _collComm;
};

#endif /*DOMAINDECOMPOSITION_H_*/
