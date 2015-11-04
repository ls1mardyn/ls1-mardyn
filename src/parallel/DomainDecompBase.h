#ifndef DOMAINDECOMPBASE_H_
#define DOMAINDECOMPBASE_H_

#include "parallel/CollectiveCommBase.h"
#include <string>

class Molecule;
class Component;
class Domain;
class ParticleContainer;
class XMLfileUnits;

typedef ParticleContainer TMoleculeContainer;

//! @brief handle boundary region and multiple processes
//! @author Martin Buchholz, Nikola Tchipev
//!
//! This program is designed to run on a HPC (High Performance Computer).
//! But sometimes one might want to execute it on a single processor, possibly even
//! without having MPI installed on that machine. One way to allow this would be to
//! have two different versions of the program, one sequential version and one parallel version.
//! But that isn't feasible, as it is hardly possible to keep them both up to date without
//! investing a lot of additional time.
//! Before describing how this problem is solved, you'll have to know a little bit about
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
//! doing all the necessary stuff for the parallelization. Further assume that
//! you have a second class which also implements this interface which is only capable
//! of handling one process but doesn't need any MPI (with only one process, there is
//! no need for message passing between processes). So the main program (or in this
//! case the class Simulation) can decide which implementation to use. When MPI is
//! available, the parallel version is used, otherwise the sequential version
class DomainDecompBase {
public:
	//! @brief The Constructor determines the own rank and the number of the neighbours                                                       */
	DomainDecompBase();

	//! @brief The Destructor finalizes MPI
	virtual ~DomainDecompBase();

	virtual void readXML(XMLfileUnits& xmlconfig);

	//! @brief exchange molecules between processes
	//!
	//! molecules which aren't in the domain of their process any
	//! more are transferred to their neighbours. Additionally, the
	//! molecules for the halo-region are transferred.
	//! This implementation is the one used in sequential mode.
	//! @param moleculeContainer needed to get those molecules which have to be exchanged
	//! @param components when creating a new Molecule-object (from the received data),
	//!                   the Molecule-constructor needs this component vector
	//! @param domain is e.g. needed to get the size of the local domain
	virtual void exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain);

	//! @brief balance the load (and optimize communication) and exchange boundary particles
	//!
	//! This method is used to perform a new decomposition of the global domain with the
	//! goal of getting the equal load (and possibly additional goals) on each process.
	//! This method is then also responsible for redistributing all particles, so after the
	//! method was called, each process has a domain with all particles belonging to this
	//! domain (as if exchangeParticles was called after the new decomposition).
	//! @param balance if true, a rebalancing should be performed, otherwise only exchange
	//! @param moleculeContainer needed for calculating load and to get the particles
	//! @param components when creating a new Molecule-object (from the received data),
	//!                   the Molecule-constructor needs this component vector
	//! @param domain is e.g. needed to get the size of the local domain
	virtual void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, Domain* domain);

	//! @brief find out whether the given position belongs to the domain of this process
	//!
	//! This method is e.g. used by a particle generator which creates particles within
	//! the bounding box (cuboid) of a local domain, but only those particles which really
	//! belong to the domain (which might not be cuboid) are allowed to be really used.
	//! @param x x-coordinate of the position to be checked
	//! @param y y-coordinate of the position to be checked
	//! @param z z-coordinate of the position to be checked
	//! @param domain might be needed to get the bounding box
	virtual bool procOwnsPos(double x, double y, double z, Domain* domain);

	void getBoundingBoxMinMax(Domain* domain, double* min, double* max);

	//! @brief get the minimum of the bounding box of this process' domain in the given dimension (0,1,2)
	//! @param dimension coordinate direction for which the minimum of the bounding box is returned
	//! @param domain here the bounding box is stored
	virtual double getBoundingBoxMin(int dimension, Domain* domain);

	//! @brief get the maximum of the bounding box of this process' domain in the given dimension (0,1,2)
	//! @param dimension coordinate direction for which the maximum of the bounding box is returned
	//! @param domain here the bounding box is stored
	virtual double getBoundingBoxMax(int dimension, Domain* domain);

	//! @brief writes information about the current decomposition into the given file
	//!        The format is not strictly defined and depends on the decomposition
	//! @param filename name of the file into which the data will be written
	//! @param domain e.g. needed to get the bounding boxes
	virtual void printDecomp(std::string filename, Domain* domain);


	//! @brief returns the own rank
	//! @return rank of the process
	virtual int getRank();

	//! @brief returns the number of processes
	//! @return number of processes
	virtual int getNumProcs();

	//! @brief synchronizes all processes
	virtual void barrier();

	//! @brief returns the time in seconds since some time in the past
	virtual double getTime();

	//! @brief returns total number of molecules
	virtual unsigned Ndistribution(unsigned localN, float* minrnd, float* maxrnd);

	//! @brief checks identity of random number generators
	virtual void assertIntIdentity(int IX);
	virtual void assertDisjunctivity(TMoleculeContainer* mm);

	//! @brief appends molecule data to the file. The format is the same as that of the input file
	//! @param filename name of the file into which the data will be written
	//! @param moleculeContainer all Particles from this container will be written to the file
	//!
	//! Currently, parallel IO isn't used.
	//! To ensure that not more than one process writes to the file at any time,
	//! there is a loop over all processes with a barrier in between
	//! @param filename name of the file into which the data will be written
	//! @param moleculeContainer all Particles from this container will be written to the file
	void writeMoleculesToFile(std::string filename, ParticleContainer* moleculeContainer);


	//##################################################################
	// The following methods with prefix "collComm" are all used
	// in the context of collective communication. Each of the methods
	// basically has to call the corresponding method from the class
	// CollectiveCommunication (or CollectiveCommDummy in the sequential
	// case). To get information about how to use this methods, read
	// the documentation of the class CollectiveCommunication.
	//##################################################################
	//! has to call init method of a CollComm class
	virtual void collCommInit(int numValues);
	//! has to call finalize method of a CollComm class
	virtual void collCommFinalize();
	//! has to call appendInt method of a CollComm class
	virtual void collCommAppendInt(int intValue);
	//! has to call appendUnsLong method of a CollComm class
	virtual void collCommAppendUnsLong(unsigned long unsLongValue);
	//! has to call appendFloat method of a CollComm class
	virtual void collCommAppendFloat(float floatValue);
	//! has to call appendDouble method of a CollComm class
	virtual void collCommAppendDouble(double doubleValue);
	//! has to call appendLongDouble method of a CollComm class
	virtual void collCommAppendLongDouble(long double longDoubleValue);
	//! has to call getInt method of a CollComm class
	virtual int collCommGetInt();
	//! has to call getUnsLong method of a CollComm class
	virtual unsigned long collCommGetUnsLong();
	//! has to call getFloat method of a CollComm class
	virtual float collCommGetFloat();
	//! has to call getDouble method of a CollComm class
	virtual double collCommGetDouble();
	//! has to call getLongDouble method of a CollComm class
	virtual long double collCommGetLongDouble();
	//! has to call allreduceSum method of a CollComm class (none in sequential version)
	virtual void collCommAllreduceSum();
	//! has to call broadcast method of a CollComm class (none in sequential version)
	virtual void collCommBroadcast(int root = 0);

private:
	virtual void handleDomainLeavingParticles(unsigned dim, ParticleContainer* moleculeContainer) const;

	virtual void populateHaloLayerWithCopies(unsigned dim, ParticleContainer* moleculeContainer) const;
private:
	CollectiveCommBase _collCommBase;
};

#endif /* DOMAINDECOMPBASE_H_ */
