/*
 * CollectiveCommunicationNonBlocking.h
 *
 *  Created on: May 2, 2016
 *      Author: seckler
 */

#pragma once

#include "CollectiveCommunicationSingleNonBlocking.h"
#include <map>
#include <utils/Logger.h>

using Log::global_log;
/**
 * CollectiveCommunicationNonBlocking provides an interface to access multiple CollectiveCommunicationSingleNonBlocking objects.
 * This allows the use of multiple different collective calls, which is needed for the different ensembles.
 * @author Steffen Seckler
 */
class CollectiveCommunicationNonBlocking {
public:
	//! Constructor, does nothing yet
	CollectiveCommunicationNonBlocking():_comms() {
	}

	//! @brief allocate memory for the values to be sent, initialize counters for specific
	//! @param key The key of the collective communication
	//! @param communicator MPI communicator for the
	//! @param numValues number of values that shall be communicated
	void init(int key, MPI_Comm communicator, int numValues) {
		if (!_comms.emplace(key).second) {
			global_log->error() << "CollectiveCommunicationNonBlocking: key "
					<< key << " already existent" << std::endl;
		}
		_comms[key].init(communicator, numValues);
	}

	//! @brief delete memory and MPI_Type
	//! @param key The key of the collective communication
	void finalize(int key) {
		_comms[key].finalize();
		_comms.erase(key);
	}

	//! Append an int value to the list of values to be sent
	//! @param key The key of the collective communication
	//! @param intValue value to be appended
	void appendInt(int key, int intValue) {
		_comms[key].appendInt(intValue);
	}

	//! Append a unsigned long value to the list of values to be sent
	//! @param key The key of the collective communication
	//! @param unsLongValue value to be appended
	void appendUnsLong(int key, unsigned long unsLongValue) {
		_comms[key].appendUnsLong(unsLongValue);
	}

	//! Append a float value to the list of values to be sent
	//! @param key The key of the collective communication
	//! @param unsLongValue value to be appended
	void appendFloat(int key, float floatValue) {
		_comms[key].appendFloat(floatValue);
	}

	//! Append a double value to the list of values to be sent
	//! @param key The key of the collective communication
	//! @param unsLongValue value to be appended
	void appendDouble(int key, double doubleValue) {
		_comms[key].appendDouble(doubleValue);
	}

	//! Append a long double value to the list of values to be sent
	//! @param key The key of the collective communication
	//! @param unsLongValue value to be appended
	void appendLongDouble(int key, long double longDoubleValue) {
		_comms[key].appendLongDouble(longDoubleValue);
	}

	//! Get the MPI communicator
	//! @param key The key of the collective communication
	//! @return MPI communicator
	MPI_Comm getTopology(int key) {
		_comms[key].getTopology();
	}

	//! Get the next value from the list, which must be int
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @param key The key of the collective communication
	//! @return the value
	int getInt(int key) {
		_comms[key].getInt();
	}

	//! Get the next value from the list, which must be unsigned long
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @param key The key of the collective communication
	//! @return the value
	unsigned long getUnsLong(int key) {
		_comms[key].getUnsLong();
	}

	//! Get the next value from the list, which must be float
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @param key The key of the collective communication
	//! @return the value
	float getFloat(int key) {
		_comms[key].getFloat();
	}

	//! Get the next value from the list, which must be double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @param key The key of the collective communication
	//! @return the value
	double getDouble(int key) {
		_comms[key].getDouble();
	}

	//! Get the next value from the list, which must be long double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @param key The key of the collective communication
	//! @return the value
	long double getLongDouble(int key) {
		_comms[key].getLongDouble();
	}

	//! Broadcast all values from the process with rank 0 to all others
	//! @param key The key of the collective communication
	//! @param root The root of the broadcast
	void broadcast(int key, int root = 0) {
		_comms[key].broadcast(root);
	}

	//! Do Allreduce off all values with reduce operation add
	//! @param key The key of the collective communication
	void allreduceSum(int key) {
		_comms[key].allreduceSum();
	}

	/**
	 * Waits until collective operation is finished.
	 * @param key The key of the collective communication
	 */
	void wait(int key) {
		_comms[key].wait();
	}

private:

	std::map<int, CollectiveCommunicationSingleNonBlocking> _comms;
};

