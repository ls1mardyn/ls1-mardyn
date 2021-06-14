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
	 */
	CollectiveCommunicationSingleNonBlocking() = default;

	/**
	 * Copy constructor is deleted, as copying the vector _tempValues will change memory addresses used in the mpi
	 * communication!
	 */
	CollectiveCommunicationSingleNonBlocking(const CollectiveCommunicationSingleNonBlocking&) = delete;

	/**
	 * Copy assignment operator is deleted, as copying the vector _tempValues will change memory addresses used in the
	 * mpi communication!
	 */
	CollectiveCommunicationSingleNonBlocking& operator=(const CollectiveCommunicationSingleNonBlocking&) = delete;

	/**
	 * Destructor
	 */
	~CollectiveCommunicationSingleNonBlocking() override {
		if (_request) {
			if (_communicationInitiated) {
				MPI_Wait(_request.get(), MPI_STATUS_IGNORE);
			}
		}
		if (_agglomeratedTypeAddOperator != MPI_OP_NULL) {
			MPI_CHECK(MPI_Op_free(&_agglomeratedTypeAddOperator));
		}
		if (_agglomeratedType != MPI_DATATYPE_NULL) {
			MPI_CHECK(MPI_Type_free(&_agglomeratedType));
		}
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
	void init(MPI_Comm communicator, int numValues, int key = 0) override {
		CollectiveCommunication::init(communicator, numValues);
#ifndef NDEBUG
		if (!_firstComm) {
			mardyn_assert(static_cast<int>(_values.size()) == numValues);
		}
#endif
	}

	// documentation in base class
	void finalize() override {
		// We intentionally use CollectiveCommBase::finalize(), as _agglomeratedType might not be MPI_DATATYPE_NULL.
		CollectiveCommBase::finalize();
		_types.clear();
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
		MPI_CHECK(MPI_Wait(_request.get(), MPI_STATUS_IGNORE));
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
		if (_agglomeratedType == MPI_DATATYPE_NULL) {
			setMPIType();
		}
		const int commutative = 1;
		valType * startOfValues = &(_tempValues[0]);
		if (_agglomeratedTypeAddOperator == MPI_OP_NULL) {
			MPI_CHECK(
					MPI_Op_create((MPI_User_function * ) CollectiveCommunication::add, commutative,
							&_agglomeratedTypeAddOperator));
		}
		MPI_CHECK(
				MPI_Iallreduce(MPI_IN_PLACE, startOfValues, 1, _agglomeratedType, _agglomeratedTypeAddOperator, _communicator, _request.get()));
		_communicationInitiated = true;
#else
		for( unsigned int i = 0; i < _types.size(); i++ ) {
			MPI_CHECK( MPI_Allreduce( MPI_IN_PLACE, &(_values[i]), 1, _types[i], MPI_SUM, _communicator ) );
		}
#endif

	}

	std::unique_ptr<MPI_Request> _request{new MPI_Request()};
	MPI_Op _agglomeratedTypeAddOperator{MPI_OP_NULL};
	bool _communicationInitiated{false};
	bool _valuesValid{false};
	/// tempValues is used for overlapped communications!
	std::vector<valType> _tempValues;
	bool _firstComm{true};
};

#endif // MPI_VERSION >= 3
