/*
 * CollectiveCommunicationSingleNonBlocking.h
 *
 *  Created on: May 2, 2016
 *      Author: seckler
 */

#pragma once

#include "CollectiveCommunication.h"
#if MPI_VERSION >= 3

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
	CollectiveCommunicationSingleNonBlocking(int key=0) {
		_key = key;
		_request = MPI_REQUEST_NULL;
		_communicationInitiated = false;
		_valuesValid = false;
		_firstComm = true;
	}

	/**
	 * Destructor
	 */
	virtual ~CollectiveCommunicationSingleNonBlocking() {
		mardyn_assert(_agglomeratedType == MPI_DATATYPE_NULL);
	}


	void allreduceSumAllowPrevious() override {
		if(!_valuesValid){
			if(_communicationInitiated){
				// if the values are not valid and communication is already initiated, first wait for the communication to finish.
				waitAndUpdateData();
			}
			// if the values were invalid, we have to do a proper allreduce!
			allreduceSum();
			_valuesValid = true;
		}

		// initiate the next allreduce
		initAllreduceSum();

	}

	// documentation in base class
	virtual void init(MPI_Comm communicator, int numValues, int key = 0) override {
		CollectiveCommBase::init(numValues);

		_communicator = communicator;
		if (_firstComm) {
			_types.reserve(numValues);
		}
#ifndef NDEBUG
		else {
			mardyn_assert(static_cast<int>(_values.size()) == numValues);
		}
#endif
	}

	// documentation in base class
	virtual void finalize() override {
		_types.clear();
		mardyn_assert(_agglomeratedType == MPI_DATATYPE_NULL);
	}

private:

	/*int testReceived() {
		int flag;
		MPI_CHECK(MPI_Test(&_request, &flag, MPI_STATUS_IGNORE));
		return flag;
	}*/

	/**
	 * Waits for communication to end and also sets the data in the correct place (values)
	 */
	void waitAndUpdateData() {
		MPI_CHECK(MPI_Wait(&_request, MPI_STATUS_IGNORE));
		// copy the temporary values to the real values!
		_values = _tempValues;

		_communicationInitiated = false;
		_valuesValid = true;
	}

	void initAllreduceSum() {
		// copy the values to a temporary vector:
		// all nonblocking operations have to be performed on this vector!
		// this is necessary to maintain the validity of the data from previous steps
		_tempValues = _values;
#if ENABLE_AGGLOMERATED_REDUCE
		setMPIType();
		MPI_Op agglomeratedTypeAddOperator;
		const int commutative = 1;
		valType * startOfValues = &(_tempValues[0]);
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
		_communicationInitiated = true;
	}

	int _key;
	MPI_Request _request;
	bool _communicationInitiated;
	bool _valuesValid;
	std::vector<valType> _tempValues;
	bool _firstComm;
};

#endif // MPI_VERSION >= 3
