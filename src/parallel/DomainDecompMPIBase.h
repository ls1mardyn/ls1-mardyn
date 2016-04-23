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

#define DIM 3

#define LOWER  0
#define HIGHER 1

class DomainDecompMPIBase: public DomainDecompBase {
public:
	DomainDecompMPIBase();
	virtual ~DomainDecompMPIBase();

	// documentation see father class (DomainDecompBase.h)
	void barrier() { MPI_CHECK( MPI_Barrier(_comm) ); }

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
	void exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType, bool removeRecvDuplicates = false);

	virtual std::vector<int> getNeighbourRanks() =0;
#if defined(ENABLE_MPI)
	virtual MPI_Comm getCommunicator(){
		return _comm;
	}
#endif
protected:

	std::vector<CommunicationPartner> _neighbours[DIM];

	MPI_Datatype _mpiParticleType;

	MPI_Comm _comm;

	//! flag, which tells whether a processor covers the whole domain along a dimension
	//! if true, we will use the methods provided by the base class for handling the
	//! respective dimension, instead of packing and unpacking messages to self
	bool _coversWholeDomain[DIM];

private:
	CollectiveCommunication _collCommunication;
};

#endif /* DOMAINDECOMPMPIBASE_H_ */
