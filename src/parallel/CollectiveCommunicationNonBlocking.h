/*
 * CollectiveCommunicationNonBlocking.h
 *
 *  Created on: May 2, 2016
 *      Author: seckler
 */
#pragma once

#include "CollectiveCommBaseInterface.h"
#include "CollectiveCommunicationSingleNonBlocking.h"
#include <map>
#include <utils/Logger.h>
#include "CollectiveCommunicationInterface.h"
#include "utils/mardyn_assert.h"

#if MPI_VERSION >= 3

/**
 * CollectiveCommunicationNonBlocking provides an interface to access multiple CollectiveCommunicationSingleNonBlocking objects.
 * This allows the use of multiple different collective calls, which is needed for the different ensembles.
 * @author Steffen Seckler
 * @note requires MPI >= 3!
 */
class CollectiveCommunicationNonBlocking: public CollectiveCommunicationInterface {
public:
	//! Constructor, does nothing yet
	CollectiveCommunicationNonBlocking() :
			_currentKey(-1), _comms() {
	}

	/**
	 * Destructor
	 */
	~CollectiveCommunicationNonBlocking() override = default;

	//! @brief allocate memory for the values to be sent, initialize counters for specific
	//! @param key The key of the collective communication, this key will be used
	//! @param communicator MPI communicator for the
	//! @param numValues number of values that shall be communicated
	void init(MPI_Comm communicator, int numValues, int key = 0) override {
		if (_currentKey != -1) {
			std::ostringstream error_message;
			error_message << "CollectiveCommunicationNonBlocking: previous communication with key " << _currentKey
					<< " not yet finalized" << std::endl;
			MARDYN_EXIT(error_message);
		}

		_currentKey = key;

		// add the key, if it is not yet existent:
		if (_comms.count(_currentKey) == 1) {
			// this happens, if the key is already existent.
			Log::global_log->debug() << "CollectiveCommunicationNonBlocking: key " << _currentKey
					<< " already existent. Reusing information." << std::endl;
		} else {
			Log::global_log->debug() << "CollectiveCommunicationNonBlocking: key " << _currentKey
								<< " not existent. Cannot reuse information." << std::endl;
			// Creates the CollectiveCommunicationSingleNonBlocking object
			auto [_, inserted] = _comms.try_emplace(_currentKey);
			if (not inserted) {
				std::ostringstream error_message;
				error_message << "CollectiveCommunicationNonBlocking: key " << _currentKey
									<< " could not be inserted. Aborting!" << std::endl;
				MARDYN_EXIT(error_message);
			}
		}
		_comms.at(_currentKey).init(communicator, numValues, _currentKey);
	}

	//! @brief delete memory and MPI_Type
	//! @param key The key of the collective communication
	void finalize() override {
		_comms.at(_currentKey).finalize();
		if (_currentKey == 0) {
			Log::global_log->debug() << "CollectiveCommunicationNonBlocking: finalizing with key " << _currentKey
					<< ", thus the entry is removed." << std::endl;
			_comms.erase(_currentKey);
		}
		_currentKey = -1;
	}

	//! Append an int value to the list of values to be sent
	//! @param intValue value to be appended
	void appendInt(int intValue) override {
		_comms.at(_currentKey).appendInt(intValue);
	}

	//! Append a unsigned long value to the list of values to be sent
	//! @param unsLongValue value to be appended
	void appendUnsLong(unsigned long unsLongValue) override {
		_comms.at(_currentKey).appendUnsLong(unsLongValue);
	}

	//! Append a float value to the list of values to be sent
	//! @param unsLongValue value to be appended
	void appendFloat(float floatValue) override {
		_comms.at(_currentKey).appendFloat(floatValue);
	}

	//! Append a double value to the list of values to be sent
	//! @param unsLongValue value to be appended
	void appendDouble(double doubleValue) override {
		_comms.at(_currentKey).appendDouble(doubleValue);
	}

	//! Append a long double value to the list of values to be sent
	//! @param unsLongValue value to be appended
	void appendLongDouble(long double longDoubleValue) override {
		_comms.at(_currentKey).appendLongDouble(longDoubleValue);
	}

	//! Get the MPI communicator
	//! @return MPI communicator
	MPI_Comm getTopology() override {
		return _comms.at(_currentKey).getTopology();
	}

	//! Get the next value from the list, which must be int
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	int getInt() override {
		return _comms.at(_currentKey).getInt();
	}

	//! Get the next value from the list, which must be unsigned long
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	unsigned long getUnsLong() override {
		return _comms.at(_currentKey).getUnsLong();
	}

	//! Get the next value from the list, which must be float
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	float getFloat() override {
		return _comms.at(_currentKey).getFloat();
	}

	//! Get the next value from the list, which must be double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	double getDouble() override {
		return _comms.at(_currentKey).getDouble();
	}

	//! Get the next value from the list, which must be long double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	long double getLongDouble() override {
		return _comms.at(_currentKey).getLongDouble();
	}

	//! Broadcast all values from the process with rank 0 to all others
	//! @param root The root of the broadcast
	void broadcast(int root = 0) override {
		_comms.at(_currentKey).broadcast(root);
	}

	//! Do Allreduce off all values with reduce operation add
	void allreduceSum() override {
		Log::global_log->debug() << "CollectiveCommunicationNonBlocking: normal Allreduce" << std::endl;
		_comms.at(_currentKey).allreduceSum();
	}

	//! Performs an all-reduce (sum), however values of previous iterations are permitted.
	//! By allowing values from previous iterations, overlapping communication is possible.
	//! One possible use case for this function is the reduction of slowly changing variables, e.g. the temperature.
	virtual void allreduceSumAllowPrevious() override {
		Log::global_log->debug() << "CollectiveCommunicationNonBlocking: nonblocking Allreduce with id "<< _currentKey << std::endl;
		mardyn_assert(_currentKey > 0); // _currentKey has to be positive non-zero and should be unique for allreduceSumAllowPrevious
		_comms.at(_currentKey).allreduceSumAllowPrevious();
	}

	void allreduceCustom(ReduceType type) override{
		_comms.at(_currentKey).allreduceCustom(type);
	}

	void scanSum() override {
		_comms.at(_currentKey).scanSum();
	}

	/**
	 * Waits until collective operation is finished.
	 * @param key The key of the collective communication
	 */
	/*
	 void waitSpecific(int key) {
	 _comms[key].waitAndUpdateData();
	 }*/

	size_t getTotalSize() override {
		size_t tmp = 0;
		for (auto& comm : _comms){
			tmp += comm.second.getTotalSize();
		}
		return tmp;
	}
private:
	int _currentKey;
	std::map<int, CollectiveCommunicationSingleNonBlocking> _comms;
};

#endif // MPI_VERSION >= 3
