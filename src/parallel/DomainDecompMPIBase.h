/*
 * DomainDecompBaseMPI.h
 *
 *  Created on: Nov 15, 2015
 *      Author: tchipevn
 */

#ifndef DOMAINDECOMPMPIBASE_H_
#define DOMAINDECOMPMPIBASE_H_

#include "DomainDecompBase.h"
#include "parallel/CollectiveCommunication.h"
#include "parallel/ParticleData.h"
#include "parallel/CommunicationPartner.h"

#include <mpi.h>
#include <vector>

#define LOWER  0
#define HIGHER 1

#define DIMgeom 3

class NeighbourCommunicationScheme;

class DomainDecompMPIBase: public DomainDecompBase {
public:
	DomainDecompMPIBase();
	virtual ~DomainDecompMPIBase();

	// documentation see father class (DomainDecompBase.h)
	void barrier() const override {
		MPI_CHECK(MPI_Barrier(_comm));
	}

	//! @brief returns total number of molecules
	unsigned Ndistribution(unsigned localN, float* minrnd, float* maxrnd);

	//! @brief checks identity of random number generators
	void assertIntIdentity(int IX);
	void assertDisjunctivity(TMoleculeContainer* mm) const override;

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
		_collCommunication.init(_comm, numValues);
	}

	void collCommFinalize() {
		_collCommunication.finalize();
	}

	void collCommAppendInt(int intValue) {
		_collCommunication.appendInt(intValue);
	}

	void collCommAppendUnsLong(unsigned long unsLongValue) {
		_collCommunication.appendUnsLong(unsLongValue);
	}

	void collCommAppendFloat(float floatValue) {
		_collCommunication.appendFloat(floatValue);
	}

	void collCommAppendDouble(double doubleValue) {
		_collCommunication.appendDouble(doubleValue);
	}

	void collCommAppendLongDouble(long double longDoubleValue) {
		_collCommunication.appendLongDouble(longDoubleValue);
	}

	int collCommGetInt() {
		return _collCommunication.getInt();
	}

	unsigned long collCommGetUnsLong() {
		return _collCommunication.getUnsLong();
	}

	float collCommGetFloat() {
		return _collCommunication.getFloat();
	}

	double collCommGetDouble() {
		return _collCommunication.getDouble();
	}

	long double collCommGetLongDouble() {
		return _collCommunication.getLongDouble();
	}

	void collCommAllreduceSum() {
		_collCommunication.allreduceSum();
	}

	void collCommBroadcast(int root = 0) {
		_collCommunication.broadcast(root);
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
			bool removeRecvDuplicates = false);

	virtual std::vector<int> getNeighbourRanks() =0;
	virtual std::vector<int> getNeighbourRanksFullShell() =0;

#if defined(ENABLE_MPI)
	MPI_Datatype getMPIParticleType() {
		return _mpiParticleType;
	}
	virtual MPI_Comm getCommunicator() {
		return _comm;
	}
#endif
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

	MPI_Comm _comm;

	NeighbourCommunicationScheme* _neighbourCommunicationScheme;
private:
	CollectiveCommunication _collCommunication;
};

#endif /* DOMAINDECOMPMPIBASE_H_ */
