#ifndef DOMAINDECOMPBASE_H_
#define DOMAINDECOMPBASE_H_

#include "parallel/CollectiveCommBase.h"
#include <string>

#ifdef ENABLE_MPI
#include "mpi.h"
#endif
#include <iostream>

#include "molecules/MoleculeForwardDeclaration.h"
#include "utils/Logger.h" // is this used?
#include "io/MemoryProfiler.h"
#include "HaloRegion.h"

class Component;
class Domain;
class ParticleContainer;
class XMLfileUnits;

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
class DomainDecompBase: public MemoryProfilable {
	friend class NeighbourCommunicationScheme;
	friend class IndirectNeighbourCommunicationScheme;
	friend class DirectNeighbourCommunicationScheme;
public:
	//! @brief The Constructor determines the own rank and the number of the neighbours                                                       */
	DomainDecompBase();

	//! @brief The Destructor finalizes MPI
	~DomainDecompBase() override;

	virtual void readXML(XMLfileUnits& xmlconfig);

	//! @brief exchange molecules between processes
	//!
	//! molecules which aren't in the domain of their process any
	//! more are transferred to their neighbours. Additionally, the
	//! molecules for the halo-region are transferred.
	//! This implementation is the one used in sequential mode.
	//! @param moleculeContainer needed to get those molecules which have to be exchanged
	//! @param domain is e.g. needed to get the size of the local domain
	void exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain);

	/**
	 * @brief Exchanges forces at the domain boundaries if it's required by the cell container.
	 * @param moleculeContainer The particle container
	 * @param domain
	 */
	virtual void exchangeForces(ParticleContainer* moleculeContainer, Domain* domain);

	/**
	 * Specifies the amount of non-blocking stages, when performing overlapping balanceAndExchange and computation.
	 * For a communication scheme, where only direct neighbours communicate, 3 stages of communication are necessary,
	 * since the particles have to be transmitted in the x-direction first, then in the y-direction, then in the z-direction.
	 * @return The amount of communication stages. Returns -1 if it is not possible.
	 */
	virtual int getNonBlockingStageCount();

	//! @brief Checks whether the balance and exchange step can be performed non-blocking.
	//!
	//! A non-blocking behaviour is typically possible, as long as no rebalancing has to be done.
	//!
	//! @param forceRebalancing if true, a rebalancing is forced;
	//! 					otherwise automatic balancing of Decomposition is applied
	//! @param moleculeContainer needed for calculating load and to get the particles
	//! @param domain is e.g. needed to get the size of the local domain
	virtual bool queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain, double etime);

	//! @brief balance the load (and optimize communication) and exchange boundary particles
	//!
	//! This method is used to perform a new decomposition of the global domain with the
	//! goal of getting the equal load (and possibly additional goals) on each process.
	//! This method is then also responsible for redistributing all particles, so after the
	//! method was called, each process has a domain with all particles belonging to this
	//! domain (as if exchangeParticles was called after the new decomposition).
	//! @param balance if true, a rebalancing is forced;
	//! 					otherwise automatic balancing of Decomposition is applied
	//! @param moleculeContainer needed for calculating load and to get the particles
	//! @param domain is e.g. needed to get the size of the local domain
	virtual void balanceAndExchange(double lastTraversalTime, bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain);

	//! @brief find out whether the given position belongs to the domain of this process
	//!
	//! This method is e.g. used by a particle generator which creates particles within
	//! the bounding box (cuboid) of a local domain, but only those particles which really
	//! belong to the domain (which might not be cuboid) are allowed to be really used.
	//! @param x x-coordinate of the position to be checked
	//! @param y y-coordinate of the position to be checked
	//! @param z z-coordinate of the position to be checked
	//! @param domain might be needed to get the bounding box
	virtual bool procOwnsPos(double x, double y, double z, Domain* domain) final;

	//! @brief get the minimum and maximum coordinate of the bounding box of this process' domain
	//! @param domain
	//! @param min lower coordinate of the bounding box
	//! @param max upper coordinate of the bounding box
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
	virtual void printDecomp(const std::string& filename, Domain* domain);


	//! @brief returns the own rank
	//! @return rank of the process
	virtual int getRank() const;

	//! @brief returns the number of processes
	//! @return number of processes
	virtual int getNumProcs() const;

	//! @brief synchronizes all processes
	virtual void barrier() const;

	//! @brief returns the time in seconds since some time in the past
	virtual double getTime() const;

	//! @brief returns total number of molecules
	virtual unsigned Ndistribution(unsigned localN, float* minrnd, float* maxrnd);

	//! @brief checks identity of random number generators
	virtual void assertIntIdentity(int IX);
	virtual void assertDisjunctivity(ParticleContainer* moleculeContainer) const;

	//! @brief returns an cutoff radius for a dimension for a global linked cells datastructure
	//!
	//! This method is e.g. used for the parallelCheckpointWriter, which builds a new global
	//! cell structure. This method returns a cutoff radius, so that each cell is fully
	//! contained in one process
	//! @param dim the dimension, which will be returned
	//! @param domain e.g. needed to get the bounding boxes
	//! @param moleculeContainer e.g. needed for the cutoff radius
	double getIOCutoffRadius(int dim, Domain* domain, ParticleContainer* moleculeContainer);


#ifdef ENABLE_MPI
	//! @brief appends molecule data to the file. The format is the same as that of the input file
	//! This version uses, MPI IO.
	//! @param filename name of the file into which the data will be written
	//! @param moleculeContainer all Particles from this container will be written to the file
	void writeMoleculesToMPIFileBinary(const std::string& filename, ParticleContainer* moleculeContainer) const;
#endif // ENABLE_MPI

	//! @brief appends molecule data to the file. The format is the same as that of the input file
	//! If MPI is enabled and binary files are supposed to be written this function will call writeMoleculesToMPIFileBinary().
	//! Otherwise, to ensure that not more than one process writes to the file at any time,
	//! there is a loop over all processes with a barrier in between.
	//! @param filename name of the file into which the data will be written
	//! @param moleculeContainer all Particles from this container will be written to the file
	//! @param binary flag, that is true if the output shall be binary
	void writeMoleculesToFile(const std::string& filename, ParticleContainer* moleculeContainer, bool binary = false) const;


	void updateSendLeavingWithCopies(bool sendTogether){
		using Log::global_log;
		// Count all processes that need to send separately
		collCommInit(1);
		collCommAppendInt(!sendTogether);
		collCommAllreduceSum();
		_sendLeavingAndCopiesSeparately = collCommGetInt();
		collCommFinalize();


		global_log->info() << "Sending leaving particles and halo copies "
				<< (sendLeavingWithCopies() ? "together" : "separately") << std::endl;
	}

	bool sendLeavingWithCopies() const{
		// No process needs to send separately => send together
		return _sendLeavingAndCopiesSeparately == 0;
	}


	//##################################################################
	// The following methods with prefix "collComm" are all used
	// in the context of collective communication. Each of the methods
	// basically has to call the corresponding method from the class
	// CollectiveCommunication (or CollectiveCommDummy in the sequential
	// case). To get information about how to use this methods, read
	// the documentation of the class CollectiveCommunication.
	//##################################################################
	//! has to call init method of a CollComm class
	virtual void collCommInit(int numValues, int key=0);
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
	//! has to call allreduceSum method of a CollComm class (none in sequential version), allows for values of previous iteration.
	virtual void collCommAllreduceSumAllowPrevious();
	//! has to call allreduceCustom method of a CollComm class (none in sequential version)
	virtual void collCommAllreduceCustom(ReduceType type);
	//! has to call scanSum method of a CollComm class (none in sequential version)
	virtual void collCommScanSum();
	//! has to call broadcast method of a CollComm class (none in sequential version)
	virtual void collCommBroadcast(int root = 0);
	//returns the ranks of the neighbours
	virtual std::vector<int> getNeighbourRanks(){
		return std::vector<int>(0);
	}
	//returns the ranks of the neighbours
	virtual std::vector<int> getNeighbourRanksFullShell(){
		std::cout << "Not yet implemented";
		return std::vector<int>(0);
	}

	//returns the ranks of all ranks
	virtual std::vector<std::vector<std::vector<int>>> getAllRanks(){
		std::cout << "Not yet implemented";
		return std::vector<std::vector<std::vector<int>>>(0);
	}
#if defined(ENABLE_MPI)
	virtual MPI_Comm getCommunicator(){
	    return MPI_COMM_WORLD;
	}
#endif

	virtual size_t getTotalSize() override {
		return _collCommBase.getTotalSize();
	}
	virtual void printSubInfo(int offset) override {
		return;
	}
	virtual std::string getName() override {
		return "DomainDecompBase";
	}

	virtual void printCommunicationPartners(std::string filename) const {};

	virtual double* getProcessTimerPointer(){
		return nullptr;
	}

protected:
	void addLeavingMolecules(std::vector<Molecule>&& invalidMolecules, ParticleContainer* moleculeContainer);

	/**
	 * Handles the sequential version of particles leaving the domain.
	 * Also used as a fall-back for the MPI variant if a process spans an entire dimension.
	 * @param dim
	 * @param moleculeContainer
	 */
	void handleDomainLeavingParticles(unsigned dim, ParticleContainer* moleculeContainer) const;

	/**
	 * Handles the sequential version of particles leaving the domain in a direct communication pattern.
	 * Hereby all particles will be moved directly to the specific region without intermediary copies.
	 * x, y and z closely resemble the offset of the current region.
	 * @param x -1, 0 or 1
	 * @param y -1, 0 or 1
	 * @param z -1, 0 or 1
	 * @param moleculeContainer
	 * @param invalidParticles used if moleculeContainer->isInvalidParticleReturner() is true
	 */
	void handleDomainLeavingParticlesDirect(const HaloRegion& haloRegion, ParticleContainer* moleculeContainer,
											std::vector<Molecule>& invalidParticles) const;

	/**
	 * @brief Does the force exchange for each dimension. Will be called for dim=0, 1 and 2.
	 * @param dim The dimension (0,1 or 2)
	 * @param moleculeContainer The particle container
	 */
	void handleForceExchange(unsigned dim, ParticleContainer* moleculeContainer) const;

	/**
	 * Does the force exchange for each direction.
	 * @param haloRegion
	 * @param moleculeContainer
	 */
	virtual void handleForceExchangeDirect(const HaloRegion& haloRegion, ParticleContainer* moleculeContainer) const;

	void populateHaloLayerWithCopies(unsigned dim, ParticleContainer* moleculeContainer) const;

	void populateHaloLayerWithCopiesDirect(const HaloRegion& haloRegion, ParticleContainer* moleculeContainer,
										   bool positionCheck = true) const;

	//! the id of the current process
	int _rank;

	//! total number of processes in the simulation
	int _numProcs;

private:
	CollectiveCommBase _collCommBase;
	int _sendLeavingAndCopiesSeparately = 0;
};

#endif /* DOMAINDECOMPBASE_H_ */
