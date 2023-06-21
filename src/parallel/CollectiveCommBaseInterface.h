#pragma once

#include <vector>

//! Reduce Types of the allreduceCustom operation
enum ReduceType {
	SUM, MIN, MAX
};

//! @brief This class is provides an interface for the base class of the collective communication.
//! @details Refactored to increase code-reuse
//! @author Steffen Seckler
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
class CollectiveCommBaseInterface {

public:

	//! @brief virtual destructor
	virtual ~CollectiveCommBaseInterface(){};

	//! Append an int value to the list of values to be stored
	//! @param intValue the value
	virtual void appendInt(int intValue) = 0;

	//! Append a unsigned long value to the list of values to be stored
	//! @param unsLongValue the value
	virtual void appendUnsLong(unsigned long unsLongValue) = 0;

	//! Append a float value to the list of values to be stored
	//! @param floatValue the value
	virtual void appendFloat(float floatValue) = 0;

	//! Append a double value to the list of values to be stored
	//! @param doubleValue the value
	virtual void appendDouble(double doubleValue) = 0;

	//! Append a long double value to the list of values to be stored
	//! @param longDoubleValue the value
	virtual void appendLongDouble(long double longDoubleValue) = 0;

	//! Performs a broadcast of the values to all processes in the communicator
	//! @param root of the broadcast
	virtual void broadcast(int /*root*/ = 0) = 0;

	//! Performs an all-reduce (sum)
	virtual void allreduceSum() = 0;

	//! Performs an allreduce operation with a custom reduce type.
	//! @param type the type of the operation
	virtual void allreduceCustom(ReduceType type) = 0;

	//! Performs a scan (sum)
	virtual void scanSum() = 0;

	//! Get the next value from the list, which must be int
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order3
	//! @return the value
	virtual int getInt() = 0;

	//! Get the next value from the list, which must be unsigned long
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	virtual unsigned long getUnsLong() = 0;

	//! Get the next value from the list, which must be float
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	virtual float getFloat() = 0;

	//! Get the next value from the list, which must be double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	virtual double getDouble() = 0;

	//! Get the next value from the list, which must be long double
	//! @details no check is performed that the type is correct,
	//! so we rely on the sanity of the programmer to ensure
	//! FIFO ordering w.r.t. append-get order
	//! @return the value
	virtual long double getLongDouble() = 0;

	//! @brief delete memory and MPI_Type
	//! @return the value
	virtual void finalize() = 0;

	//! get the size of the entire object allocated memory
	//! @return size of the dynamically allocated memory in bytes
	virtual size_t getTotalSize() = 0;
};

