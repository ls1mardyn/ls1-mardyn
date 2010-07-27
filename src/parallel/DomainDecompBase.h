#ifndef DOMAINDECOMPBASE_H_
#define DOMAINDECOMPBASE_H_

#include <vector>
#include <iostream>
#include <map>

class Molecule;
class Component;
class Domain;
class ParticleContainer;

typedef ParticleContainer TMoleculeContainer;

//! @brief handle boundary region and multiple processes
//! @author Martin Buchholz
//!
//! This program is designed to run on a HPC (High Performance Computer).
//! But sometimes one might want to execute it on a single processor, possibly even
//! without having MPI installed on that machine. One way to allow this would be to
//! have two different versions of the program, one sequential version and one parallel version.
//! But that isn't feasible, as it is hardly possible to keep them both up to date without
//! investing a lot of additional time.
//! Befor describing how this problem is solved, you'll have to know a little bit about
//! how the parallel version works.
//!
//! At the moment, domain decomposition is used. The basic idea is, that the region which
//! is simulated (usually a cube) is divided into several smaller regions. Each process
//! works on one of those regions. But as the molecules move between regions, the processes
//! have to communicate to exchange molecules. Also for the calculation of some values
//! (e.g. macroscopic values), communication is necessary. Each time processes need
//! to communicate, a method of this interface (that means a method of a class implementing
//! this interface) is called, which then somehow does the communication.
//!
//! Assume you have an class which implements this interface and which is capable of
//! doing all the necessary stuff for the parallelisation. Further assume that
//! you have a second class which also implements this interface which is only capable
//! of handling one process but doesn't need any MPI (with only one process, there is
//! no need for message passing between processes). So the main program (or in this
//! case the class Simulation) can decide which implementation to use. When MPI is
//! available, the parallel version is used, otherwise the sequential version
class DomainDecompBase {
public:
	//! @brief The Constructor determines the own rank and the number of the neighbours                                                       */
	DomainDecompBase() {
	}

	//! @brief The Destructor finalizes MPI
	virtual ~DomainDecompBase() {
	}

	//! @brief exchange molecules between processes
	//!
	//! molecules which aren't in the domain of their process any
	//! more are transferred to their neighbours. Additionally, the
	//! molecules for the halo-region are transferred.
	//! @param moleculeContainer needed to get those molecules which have to be exchanged
	//! @param components when creating a new Molecule-object (from the recieved data),
	//!                   the Molecule-constructor needs this component vector
	//! @param domain is e.g. needed to get the size of the local domain
	virtual void exchangeMolecules(ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain) = 0;

	//! @brief balance the load (and optimise communication) and exchange boundary particles
	//!
	//! This method is used to perform a new decomposition of the global domain with the
	//! goal of getting the equal load (and possibly additional goals) on each process.
	//! This method is then also responsible for redistributing all particles, so after the
	//! method was called, each process has a domain with all particles belonging to this
	//! domain (as if exchangeParticles was called after the new decomposition).
	//! @param balance if true, a rebalancing should be performed, otherwise only exchange
	//! @param moleculeContainer needed for calculating load and to get the particles
	//! @param components when creating a new Molecule-object (from the recieved data),
	//!                   the Molecule-constructor needs this component vector
	//! @param domain is e.g. needed to get the size of the local domain
	virtual void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain) = 0;

	//! @brief find out whether the given position belongs to the domain of this process
	//!
	//! This method is e.g. used by a particle generator which creates particles within
	//! the bounding box (cuboid) of a local domain, but only those particles which really
	//! belong to the domain (which might not be cuboid) are allowed to be really used.
	//! @param x x-coordinate of the position to be checked
	//! @param y y-coordinate of the position to be checked
	//! @param z z-coordinate of the position to be checked
	//! @param domain might be needed to get the bounding box
	virtual bool procOwnsPos(double x, double y, double z, Domain* domain) = 0;

	//! @brief returns a guaranteed distance of (x,y,z) to the local domain
	//!
	//! This method is e.g. used by a particle generator which creates clusters (nuclei, drops).
	//! Only if the cluster is close (cluster radius larger then guaranteedDistance) to the
	//! domain the particles have to be created.
	//! @param x x-coordinate of the position for which the guaranteed distance is returned
	//! @param y y-coordinate of the position for which the guaranteed distance is returned
	//! @param z z-coordinate of the position for which the guaranteed distance is returned
	//! @param domain might be needed to get the bounding box
	virtual double guaranteedDistance(double x, double y, double z, Domain* domain) = 0;

	//! @brief counts the number of molecules of each component type.
	//!
	//! This method is usually only needed once in the beginning of the simulation and only
	//! if the particles were not read in from a single file but read in from one file per proc or
	//! if the particles were created by each proc seperately.
	//! @param moleculeContainer container for the molecules
	//! @param compCount vector which has to have the size which equals the number of components
	//!                  this method will will the vector with the number of molecules for each
	//!                  of the components (in the global domain)
	//! @return the number of molecules in the global domain is returned
	virtual unsigned long countMolecules(ParticleContainer* moleculeContainer, std::vector<unsigned long> &compCount) = 0;

	//! @brief get the minimum of the bounding box of this process' domain in the given dimension (0,1,2)
	//! @param dimension coordinate direction for which the minimum of the bounding box is returned
	//! @param domain here the bounding box is stored
	virtual double getBoundingBoxMin(int dimension, Domain* domain) = 0;

	//! @brief get the maximum of the bounding box of this process' domain in the given dimension (0,1,2)
	//! @param dimension coordinate direction for which the maximum of the bounding box is returned
	//! @param domain here the bounding box is stored
	virtual double getBoundingBoxMax(int dimension, Domain* domain) = 0;

	//! @brief writes information about the current decomposition into the given file
	//!        The format is not strictly defined and depends on the decompositio
	//! @param filename name of the file into which the data will be written
	//! @param domain e.g. needed to get the bounding boxes
	virtual void printDecomp(std::string filename, Domain* domain) = 0;

	//! @brief appends molecule data to the file. The format is the same as that of the input file
	//! @param filename name of the file into which the data will be written
	//! @param moleculeContainer all Particles from this container will be written to the file
	virtual void writeMoleculesToFile(std::string filename, ParticleContainer* moleculeContainer) = 0;

	//! @brief returns the own rank
	//! @return rank of the process
	virtual int getRank() = 0;

	//! @brief returns the number of processes
	//! @return number of processes
	virtual int getNumProcs() = 0;

	//! @brief synchronises all processes
	virtual void barrier() = 0;

	//! @brief returns the time in seconds since some time in the past
	virtual double getTime() = 0;

	//! @brief returns total number of molecules
	virtual unsigned Ndistribution(unsigned localN, float* minrnd, float* maxrnd) = 0;

	//! @brief checks identity of random number generators
	virtual void assertIntIdentity(int IX) = 0;
	virtual void assertDisjunctivity(TMoleculeContainer* mm) = 0;

	//##################################################################
	// The following methods with prefix "collComm" are all used
	// in the context of collective communication. Each of the methods
	// basically has to call the corresponding method from the class
	// CollectiveCommunication (or CollectiveCommDummy in the sequential
	// case). To get information about how to use this methods, read
	// the documentation of the class CollectiveCommunication.
	//##################################################################
	//! has to call init method of a CollComm class
	virtual void collCommInit(int numValues) = 0;
	//! has to call finalize method of a CollComm class
	virtual void collCommFinalize() = 0;
	//! has to call appendInt method of a CollComm class
	virtual void collCommAppendInt(int intValue) = 0;
	//! has to call appendUnsLong method of a CollComm class
	virtual void collCommAppendUnsLong(unsigned long unsLongValue) = 0;
	//! has to call appendFloat method of a CollComm class
	virtual void collCommAppendFloat(float floatValue) = 0;
	//! has to call appendDouble method of a CollComm class
	virtual void collCommAppendDouble(double doubleValue) = 0;
	//! has to call appendLongDouble method of a CollComm class
	virtual void collCommAppendLongDouble(long double longDoubleValue) = 0;
	//! has to call getInt method of a CollComm class
	virtual int collCommGetInt() = 0;
	//! has to call getUnsLong method of a CollComm class
	virtual unsigned long collCommGetUnsLong() = 0;
	//! has to call getFloat method of a CollComm class
	virtual float collCommGetFloat() = 0;
	//! has to call getDouble method of a CollComm class
	virtual double collCommGetDouble() = 0;
	//! has to call getLongDouble method of a CollComm class
	virtual long double collCommGetLongDouble() = 0;
	//! has to call allreduceSum method of a CollComm class (none in sequential version)
	virtual void collCommAllreduceSum() = 0;
	//! has to call broadcast method of a CollComm class (none in sequential version)
	virtual void collCommBroadcast() = 0;
};

#endif /*DOMAINDECOMPBASE_H_*/
