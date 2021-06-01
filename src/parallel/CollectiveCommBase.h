#pragma once

#include <vector>
#include <stddef.h>
#include "CollectiveCommBaseInterface.h"
//! @brief This class is a dummy class which ensures that the collective communication
//!        commands also work if the program is executed sequentially without MPI.
//! @details Refactored to increase code-reuse
//! @author Nikola Tchipev, Martin Buchholz
//!
//! To fully understand the purpose of this class, read the documentation
//! of the class CollectiveCommunication. The application developer should not
//! have to care whether the program will be executed in parallel or sequentially.
//! The only thing that the application developer should know that there might
//! be several processes and at some point in the simulation those have to do a
//! collective communication to broadcast or reduce values. Those values are 
//! given to the domain decomposition method, the communication command is executed
//! by the domain decomp and the values can be read again. This class just has to
//! ensure that this also works if there is no real domain decomposition but the
//! class DomainDecompBase is used.
class CollectiveCommBase : public virtual CollectiveCommBaseInterface{

protected:
	//! As in C++ arrays have to contain only one type of variable,
	//! but this class (see documentation of CollectiveCommunication class)
	//! shall be able to store values of different types for the transfer,
	//! this union is needed to be able to store the values of different
	//! types in one array.
	union valType {
		int v_int;
		unsigned long v_unsLong;
		float v_float;
		double v_double;
		long double v_longDouble;
	};

public:
	virtual ~CollectiveCommBase() {
	}

	//! @brief allocate memory for the values to be stored, initialize getter-iterator
	//! @param numValues number of values that shall be stored
	void init(int numValues) {
		_values.reserve(numValues);
		_getter = _values.begin();
	}

	//! Append an int value to the list of values to be stored
	//! @param intValue the value
	virtual void appendInt(int intValue) override {
		valType toPush;
		toPush.v_int = intValue;
		_values.push_back(toPush);
	}

	//! Append a unsigned long value to the list of values to be stored
	//! @param unsLongValue the value
	virtual void appendUnsLong(unsigned long unsLongValue) override {
		valType toPush;
		toPush.v_unsLong = unsLongValue;
		_values.push_back(toPush);
	}

	//! Append a float value to the list of values to be stored
	//! @param floatValue the value
	virtual void appendFloat(float floatValue) override {
		valType toPush;
		toPush.v_float = floatValue;
		_values.push_back(toPush);
	}

	//! Append a double value to the list of values to be stored
	//! @param doubleValue the value
	virtual void appendDouble(double doubleValue) override {
		valType toPush;
		toPush.v_double = doubleValue;
		_values.push_back(toPush);
	}

	//! Append a long double value to the list of values to be stored
	//! @param longDoubleValue the value
	virtual void appendLongDouble(long double longDoubleValue) override {
		valType toPush;
		toPush.v_longDouble = longDoubleValue;
		_values.push_back(toPush);
	}

	//! Performs a broadcast of the values to all processes in the communicator
	//! @param root of the broadcast
	virtual void broadcast(int /*root*/ = 0) override {
	}

	//! Performs an all-reduce (sum)
	virtual void allreduceSum() override {
	}

	// doku in base class
	virtual void allreduceCustom(ReduceType type) override{
	}

	//! Performs a scan (sum)
	virtual void scanSum() override {
	}

	//! Get the next value from the list, which must be int
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order3
	//! @return the value
	virtual int getInt() override {
		return (_getter++)->v_int;
	}

	//! Get the next value from the list, which must be unsigned long
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	virtual unsigned long getUnsLong() override {
		return (_getter++)->v_unsLong;
	}

	//! Get the next value from the list, which must be float
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	virtual float getFloat() override {
		return (_getter++)->v_float;
	}

	//! Get the next value from the list, which must be double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	virtual double getDouble() override {
		return (_getter++)->v_double;
	}

	//! Get the next value from the list, which must be long double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	virtual long double getLongDouble() override {
		return (_getter++)->v_longDouble;
	}

	//! @brief delete memory and MPI_Type
	//! @return the value
	virtual void finalize() override {
		_values.clear();
	}

	virtual size_t getTotalSize() override {
		return _values.capacity() * sizeof(valType);
	}

protected:
	//! Vector to store the values which shall be communicated
	std::vector<valType> _values;

	//! Iterator to extract the values which were communicated
	std::vector<valType>::const_iterator _getter;

};

