/*
 * CollectiveCommunicationSingleNonBlocking.h
 *
 *  Created on: May 2, 2016
 *      Author: seckler
 */

#pragma once

#include "CollectiveCommunication.h"

/**
 * CollectiveCommunicationSingleNonBlocking extends the CollectiveCommunication by an implementation of nonblocking collectives.
 * @author Steffen Seckler
 */
class CollectiveCommunicationSingleNonBlocking: public CollectiveCommunication {

public:
	/**
	 * Constructor
	 * @param key the key of the collective operation
	 */
	CollectiveCommunicationSingleNonBlocking(int key) {
		_key = key;
		_request = 0;
	}

	/**
	 * Destructor
	 */
	virtual ~CollectiveCommunicationSingleNonBlocking() {
		mardyn_assert(_agglomeratedType == MPI_DATATYPE_NULL);
	}

	void initAllreduceSum(int /*key*/) {
#if ENABLE_AGGLOMERATED_REDUCE
		setMPIType();
		MPI_Op agglomeratedTypeAddOperator;
		const int commutative = 1;
		valType * startOfValues = &(_values[0]);
		MPI_CHECK(
				MPI_Op_create(
						(MPI_User_function * ) CollectiveCommunication::add,
						commutative, &agglomeratedTypeAddOperator));
		MPI_CHECK(
				MPI_Iallreduce(MPI_IN_PLACE, startOfValues, 1, _agglomeratedType, agglomeratedTypeAddOperator, _communicator, & _request));
		MPI_CHECK(MPI_Op_free(&agglomeratedTypeAddOperator));
		MPI_CHECK(MPI_Type_free(&_agglomeratedType));
#else
		for( int i = 0; i < _numValues; i++ ) {
			MPI_CHECK( MPI_Allreduce( MPI_IN_PLACE, &(_values[i]), 1, _types[i], MPI_SUM, _communicator ) );
		}
#endif
	}

	/**
	 * waits until collective operation is finished.
	 */
	void wait(){
		MPI_CHECK(MPI_Wait(&_request, MPI_STATUS_IGNORE));
		mardyn_assert(testReceived());
	}

	// documented in base class
	int getInt() {
		mardyn_assert(testReceived());
		return CollectiveCommBase::getInt();
	}

	// documented in base class
	unsigned long getUnsLong() {
		mardyn_assert(testReceived());
		return CollectiveCommBase::getInt();
	}

	// documented in base class
	float getFloat() {
		mardyn_assert(testReceived());
		return CollectiveCommBase::getInt();
	}

	// documented in base class
	double getDouble() {
		mardyn_assert(testReceived());
		return CollectiveCommBase::getInt();
	}

	// documented in base class
	long double getLongDouble() {
		mardyn_assert(testReceived());
		return CollectiveCommBase::getInt();
	}
private:
	int testReceived() {
		int flag;
		MPI_CHECK(MPI_Test(&_request, &flag, MPI_STATUS_IGNORE));
		return flag;
	}
	int _key;
	MPI_Request _request;
};
