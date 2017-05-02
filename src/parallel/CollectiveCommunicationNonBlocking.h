/*
 * CollectiveCommunicationNonBlocking.h
 *
 *  Created on: May 2, 2016
 *      Author: seckler
 */
#pragma once

#include "CollectiveCommunicationInterface.h"
#include "CollectiveCommBaseInterface.h"
#include "CollectiveCommunicationSingleNonBlocking.h"
#include <map>
#include <utils/Logger.h>

using Log::global_log;
/**
 * CollectiveCommunicationNonBlocking provides an interface to access multiple CollectiveCommunicationSingleNonBlocking objects.
 * This allows the use of multiple different collective calls, which is needed for the different ensembles.
 * @author Steffen Seckler
 */
class CollectiveCommunicationNonBlocking: public CollectiveCommunicationInterface, public CollectiveCommBaseInterface {
public:
	//! Constructor, does nothing yet
	CollectiveCommunicationNonBlocking() :
			_comms(), _currentKey(-1) {
	}
	/**
	 * Destructor
	 */
	virtual ~CollectiveCommunicationNonBlocking() {
	}

	//! @brief allocate memory for the values to be sent, initialize counters for specific
	//! @param key The key of the collective communication, this key will be used
	//! @param communicator MPI communicator for the
	//! @param numValues number of values that shall be communicated
	void init(MPI_Comm communicator, int numValues, int key = 0) override {
		if(_currentKey != -1){
			global_log->error() << "CollectiveCommunicationNonBlocking: previous communication with key " << _currentKey
					<< " not yet finalized" << std::endl;
			Simulation::exit(234);
		}
		_currentKey = key;
		if (!_comms.emplace(key).second) {
			// this happens, if the key is already existent.
			global_log->debug() << "CollectiveCommunicationNonBlocking: key "
					<< key << " already existent. Reusing information." << std::endl;
		}
		else{
			_comms[key].init(communicator, numValues);
		}
	}

	//! @brief delete memory and MPI_Type
	//! @param key The key of the collective communication
	void finalize() override {
		if(_currentKey==0){
			global_log->debug() << "CollectiveCommunicationNonBlocking: finalizing with key "
								<< _currentKey << ", thus the entry is removed." << std::endl;
			_comms[_currentKey].finalize();
			_comms.erase(_currentKey);
		}
		_currentKey = -1;
	}

	//! Append an int value to the list of values to be sent
	//! @param intValue value to be appended
	void appendInt(int intValue) override {
		_comms[_currentKey].appendInt(intValue);
	}

	//! Append a unsigned long value to the list of values to be sent
	//! @param unsLongValue value to be appended
	void appendUnsLong(unsigned long unsLongValue) override {
		_comms[_currentKey].appendUnsLong(unsLongValue);
	}

	//! Append a float value to the list of values to be sent
	//! @param unsLongValue value to be appended
	void appendFloat(float floatValue) override {
		_comms[_currentKey].appendFloat(floatValue);
	}

	//! Append a double value to the list of values to be sent
	//! @param unsLongValue value to be appended
	void appendDouble(double doubleValue) override {
		_comms[_currentKey].appendDouble(doubleValue);
	}

	//! Append a long double value to the list of values to be sent
	//! @param unsLongValue value to be appended
	void appendLongDouble(long double longDoubleValue) override {
		_comms[_currentKey].appendLongDouble(longDoubleValue);
	}

	//! Get the MPI communicator
	//! @return MPI communicator
	MPI_Comm getTopology() override {
		return _comms[_currentKey].getTopology();
	}

	//! Get the next value from the list, which must be int
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	int getInt() override {
		return _comms[_currentKey].getInt();
	}

	//! Get the next value from the list, which must be unsigned long
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	unsigned long getUnsLong() override {
		return _comms[_currentKey].getUnsLong();
	}

	//! Get the next value from the list, which must be float
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	float getFloat() override {
		return _comms[_currentKey].getFloat();
	}

	//! Get the next value from the list, which must be double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	double getDouble() override {
		return _comms[_currentKey].getDouble();
	}

	//! Get the next value from the list, which must be long double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	long double getLongDouble() override {
		return _comms[_currentKey].getLongDouble();
	}

	//! Broadcast all values from the process with rank 0 to all others
	//! @param root The root of the broadcast
	void broadcast(int root = 0) override {
		_comms[_currentKey].broadcast(root);
	}

	//! Do Allreduce off all values with reduce operation add
	void allreduceSum() override {
		_comms[_currentKey].allreduceSum();
	}

	//! Performs an all-reduce (sum), however values of previous iterations are permitted.
	//! By allowing values from previous iterations, overlapping communication is possible.
	//! One possible use case for this function is the reduction of slowly changing variables, e.g. the temperature.
	virtual void allreduceSumAllowPrevious() override {
		mardyn_assert(_currentKey > 0);  // _currentKey has to be positive non-zero and should be unique for allreduceSumAllowPrevious
		_comms[_currentKey].allreduceSumAllowPrevious();
	}

	/**
	 * Waits until collective operation is finished.
	 * @param key The key of the collective communication
	 */
	/*
	void waitSpecific(int key) {
		_comms[key].waitAndUpdateData();
	}*/

private:
	int _currentKey;
	std::map<int, CollectiveCommunicationSingleNonBlocking> _comms;
};

