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
	CollectiveCommunicationSingleNonBlocking(int key = 0) {
		_request = nullptr;
		_key = key;
		_communicationInitiated = false;
		_valuesValid = false;
		_firstComm = true;
	}

	void instantiate() {
		_request = new MPI_Request();
	}

	void destroy() {
		if (_request) {
			if (_communicationInitiated) {
				MPI_Wait(_request, MPI_STATUS_IGNORE);
			}
			delete _request;
		}
	}
	/**
	 * Destructor
	 */
	virtual ~CollectiveCommunicationSingleNonBlocking() {
		mardyn_assert(_agglomeratedType == MPI_DATATYPE_NULL);
	}

	void allreduceSumAllowPrevious() override {
		if (!_valuesValid) {
			if (_communicationInitiated) {
				// if the values are not valid and communication is already initiated, first wait for the communication to finish.
				waitAndUpdateData();
			}
			// if the values were invalid, we have to do a proper allreduce!
			allreduceSum();
			_tempValues = _values;  // save it somewhere safe for next iteration!
			_valuesValid = true;
		} else {

			if (_communicationInitiated) {
				// if the values are not valid and communication is already initiated, first wait for the communication to finish.
				// safe the data, that needs to be sent
				std::vector<valType> toSendValues = _values;
				// get the new data, this updates _values
				waitAndUpdateData();
				// initiate the next allreduce
				initAllreduceSum(toSendValues);
			} else {
				std::vector<valType> previous = _tempValues;
				// initiate the next allreduce
				initAllreduceSum(_values);  // this updates _tempValues
				// restore saved data from previous operations.
				_values = previous;
			}

		}
	}

	// documentation in base class
	virtual void init(MPI_Comm communicator, int numValues, int key = 0) override {
		CollectiveCommunication::init(communicator, numValues);
#ifndef NDEBUG
		if (!_firstComm) {
			mardyn_assert(static_cast<int>(_values.size()) == numValues);
		}
#endif
	}

	// documentation in base class
	virtual void finalize() override {
		CollectiveCommunication::finalize();
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
		MPI_CHECK(MPI_Wait(_request, MPI_STATUS_IGNORE));
		// copy the temporary values to the real values!
		_values = _tempValues;

		_communicationInitiated = false;
		_valuesValid = true;
	}

	void initAllreduceSum(std::vector<valType> values) {
		// copy the values to a temporary vector:
		// all nonblocking operations have to be performed on this vector!
		// this is necessary to maintain the validity of the data from previous steps
		_tempValues = values;
#if ENABLE_AGGLOMERATED_REDUCE
		setMPIType();
		MPI_Op agglomeratedTypeAddOperator;
		const int commutative = 1;
		valType * startOfValues = &(_tempValues[0]);
		MPI_CHECK(
				MPI_Op_create((MPI_User_function * ) CollectiveCommunication::add, commutative,
						&agglomeratedTypeAddOperator));
		MPI_CHECK(
				MPI_Iallreduce(MPI_IN_PLACE, startOfValues, 1, _agglomeratedType, agglomeratedTypeAddOperator, _communicator, _request));
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
	MPI_Request* _request;
	bool _communicationInitiated;
	bool _valuesValid;
	std::vector<valType> _tempValues;
	bool _firstComm;
};

#endif // MPI_VERSION >= 3
