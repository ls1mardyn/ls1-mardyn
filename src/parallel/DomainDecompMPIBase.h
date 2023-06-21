/*
 * DomainDecompBaseMPI.h
 *
 *  Created on: Nov 15, 2015
 *      Author: tchipevn
 */

#ifndef DOMAINDECOMPMPIBASE_H_
#define DOMAINDECOMPMPIBASE_H_

#include <mpi.h>
#include <vector>
#include <memory>

#include "utils/Logger.h"
#include "DomainDecompBase.h"
#include "CollectiveCommunicationInterface.h"
#include "CommunicationPartner.h"
#include "ParticleDataForwardDeclaration.h"

#define LOWER  0
#define HIGHER 1

#define DIMgeom 3

class NeighbourCommunicationScheme;

struct HaloRegion;

class DomainDecompMPIBase: public DomainDecompBase {
public:
	DomainDecompMPIBase();
	virtual ~DomainDecompMPIBase();

	// documentation see father class (DomainDecompBase.h)
	void barrier() const override {
		MPI_CHECK(MPI_Barrier(_comm));
	}

	//! @brief returns total number of molecules
	unsigned Ndistribution(unsigned localN, float* minrnd, float* maxrnd) override;

	//! @brief checks identity of random number generators
	void assertIntIdentity(int IX) override;
	void assertDisjunctivity(ParticleContainer* moleculeContainer) const override;

	//##################################################################
	// The following methods with prefix "collComm" are all used
	// in the context of collective communication. Each of the methods
	// basically has to call the corresponding method from the class
	// CollectiveCommunication (or CollectiveCommDummy in the sequential
	// case). To get information about how to use this methods, read
	// the documentation of the class CollectiveCommunication and of the
	// father class of this class (DomainDecompBase.h)
	//##################################################################
	void collCommInit(int numValues, int key=0) override {
		_collCommunication->init(_comm, numValues, key);
	}

	void collCommFinalize() override {
		_collCommunication->finalize();
	}

	void collCommAppendInt(int intValue) override {
		_collCommunication->appendInt(intValue);
	}

	void collCommAppendUnsLong(unsigned long unsLongValue) override {
		_collCommunication->appendUnsLong(unsLongValue);
	}

	void collCommAppendFloat(float floatValue) override {
		_collCommunication->appendFloat(floatValue);
	}

	void collCommAppendDouble(double doubleValue) override {
		_collCommunication->appendDouble(doubleValue);
	}

	void collCommAppendLongDouble(long double longDoubleValue) override {
		_collCommunication->appendLongDouble(longDoubleValue);
	}

	int collCommGetInt() override {
		return _collCommunication->getInt();
	}

	unsigned long collCommGetUnsLong() override {
		return _collCommunication->getUnsLong();
	}

	float collCommGetFloat() override {
		return _collCommunication->getFloat();
	}

	double collCommGetDouble() override {
		return _collCommunication->getDouble();
	}

	long double collCommGetLongDouble() override {
		return _collCommunication->getLongDouble();
	}

	void collCommAllreduceSum() override {
		_collCommunication->allreduceSum();
	}

	void collCommAllreduceSumAllowPrevious() override;

	void collCommAllreduceCustom(ReduceType type) override {
		_collCommunication->allreduceCustom(type);
	}

	void collCommScanSum() override {
		_collCommunication->scanSum();
	}

	void collCommBroadcast(int root = 0) override {
		_collCommunication->broadcast(root);
	}

	/**
	 * Initialises the non-blocking balance and exchange.
	 * Nothing really important needs to be done here.
	 * Some ideas: decide between possible communication schemes,...
	 *
	 * @param forceRebalancing true if rebalancing should be forced
	 * @param moleculeContainer pointer to the molecule container
	 * @param domain pointer to the domain
	 */
	virtual void balanceAndExchangeInitNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer,
			Domain* domain);

	/**
	 * Prepares the stageNumber'th stage.
	 * This includes sending particles, that have to be send in that stage.
	 * @param forceRebalancing true if rebalancing should be forced
	 * @param moleculeContainer pointer to the molecule container
	 * @param domain pointer to the domain
	 * @param stageNumber the number of the stage, the communication is in.
	 */
	virtual void prepareNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber) = 0;

	/**
	 * Finishes the stageNumber'th stage.
	 * This includes receiving the particles, that have to be received in that stage.
	 * @param forceRebalancing true if rebalancing should be forced
	 * @param moleculeContainer pointer to the molecule container
	 * @param domain pointer to the domain
	 * @param stageNumber the number of the stage, the communication is in.
	 */
	virtual void finishNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber) = 0;

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
	//! @param domain is e.g. needed to get the size of the local domain
	void exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
							  bool doHaloPositionCheck = true, bool removeRecvDuplicates = false);

	void exchangeForces(ParticleContainer* moleculeContainer, Domain* domain) override;

	std::vector<int> getNeighbourRanks() override = 0;
	std::vector<int> getNeighbourRanksFullShell() override = 0;

	virtual std::vector<CommunicationPartner> getNeighboursFromHaloRegion(Domain* domain, const HaloRegion& haloRegion, double cutoff) = 0;


#if defined(ENABLE_MPI)
	MPI_Datatype getMPIParticleType() {
		return _mpiParticleType;
	}
	MPI_Datatype getMPIParticleForceType() {
		return _mpiParticleForceType;
	}
	MPI_Comm getCommunicator() override {
		return _comm;
	}
#endif

	/** @brief Read in XML configuration for DomainDecompMPIBase.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <parallelisation type="DomainDecomposition" OR "KDDecomposition">
	     <!--default: indirect, unless in autopas mode-->
	   	 <CommunicationScheme>indirect OR direct OR direct-pp</CommunicationScheme>
	   	 <!--default: no-->
	   	 <overlappingCollectives>yes OR no</overlappingCollectives>
	   	 <!--default: 5-->
	   	 <overlappingStartAtStep></overlappingStartAtStep>
	   	 <!--default: yes-->
	   	 <useSequentialFallback>yes OR no</useSequentialFallback>
	     <!-- structure handled by DomainDecomposition or KDDecomposition -->
	   </parallelisation>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig) override;

	/**
	 * Sets the communicationScheme.
	 * @note If this function is called dynamically, make sure to reinitialise the CommunicationPartners before exchanging molecules!
	 * @param scheme
	 * @param comScheme
	 */
	virtual void setCommunicationScheme(const std::string& scheme, const std::string& comScheme);

	// documentation in base class
	virtual int getNonBlockingStageCount() override;

	virtual size_t getTotalSize() override;

	virtual void printSubInfo(int offset) override;

    virtual void printDecomp(const std::string &filename, Domain *domain, ParticleContainer *particleContainer) override;

	virtual std::string getName() override {
		return "DomainDecompMPIBase";
	}

	void printCommunicationPartners(std::string filename) const override;
protected:

	/**
	 * Prepares the stageNumber'th stage.
	 * This includes sending particles, that have to be send in that stage.
	 * @param moleculeContainer pointer to the molecule container
	 * @param domain pointer to the domain
	 * @param stageNumber the number of the stage, the communication is in.
	 */
	virtual void prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates = false);

	/**
	 * Finishes the stageNumber'th stage.
	 * This includes receiving the particles, that have to be received in that stage.
	 * @param forceRebalancing true if rebalancing should be forced
	 * @param moleculeContainer pointer to the molecule container
	 * @param domain pointer to the domain
	 * @param stageNumber the number of the stage, the communication is in.
	 */
	virtual void finishNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates = false);

	MPI_Datatype _mpiParticleType;
	MPI_Datatype _mpiParticleForceType;

	MPI_Comm _comm;

	std::unique_ptr<NeighbourCommunicationScheme> _neighbourCommunicationScheme;

	/**
	 * Indicates whether the direct-pp communication scheme should be forced.
	 * Flag to indicate necessity to use the direct-pp neighbour communication scheme. Can be overwritten by any child class.
	 * @note Intended to be used in `readXML()`.
	 */
	bool _forceDirectPP{false};
private:
	std::unique_ptr<CollectiveCommunicationInterface> _collCommunication;
	/// Defines when to start the overlapping collective communication.
	/// In the very first steps, the thermostats heavily influence the simulation.
	/// To prevent large deviations from a simulation without overlapping collectives, overlapping collectives can be disabled at the start.
	/// Typically, around five steps are reasonable for this.
	unsigned long _overlappingStartAtStep {5ul};
};

#endif /* DOMAINDECOMPMPIBASE_H_ */
