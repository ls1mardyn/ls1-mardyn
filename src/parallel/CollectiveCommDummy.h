#ifndef COLLECTIVECOMMDUMMY_H_
#define COLLECTIVECOMMDUMMY_H_


//! @brief This class is just a dummy class which ensures that the collective communication
//!        commands also work if the program is executed sequentially without MPI
//! @author Martin Buchholz
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
//! class DomainDecompDummy is used.
class CollectiveCommDummy {


	//! Read documentation for this union in CollectiveCommunication.h
	union valType {
		int val_int;
		unsigned long val_unsLong;
		float val_float;
		double val_double;
		long double val_longDouble;
	};


public:
	//! @brief allocate memory for the values to be stored, initialize counters 
	//! @param numValues number of values that shall be stored
	void init(int numValues) {
		_values = new valType[numValues];
		_numValues = numValues;
		_setCounter = 0;
		_getCounter = 0;
	}

	//! @brief delete memory and MPI_Type
	void finalize() {
		delete[] _values;
	}

	//! Append an int value to the list of values to be stored
	void appendInt(int intValue) {
		_values[_setCounter++].val_int = intValue;
	}
	//! Append a unsigned long value to the list of values to be stored
	void appendUnsLong(unsigned long unsLongValue) {
		_values[_setCounter++].val_unsLong = unsLongValue;
	}
	//! Append a float value to the list of values to be stored
	void appendFloat(float floatValue) {
		_values[_setCounter++].val_float = floatValue;
	}
	//! Append a double value to the list of values to be stored
	void appendDouble(double doubleValue) {
		_values[_setCounter++].val_double = doubleValue;
	}
	//! Append a long double value to the list of values to be stored
	void appendLongDouble(double longDoubleValue) {
		_values[_setCounter++].val_longDouble = longDoubleValue;
	}

	//! Get the next value from the list, which must be int
	int getInt() {
		return _values[_getCounter++].val_int;
	}
	//! Get the next value from the list, which must be unsigned long
	unsigned long getUnsLong() {
		return _values[_getCounter++].val_unsLong;
	}
	//! Get the next value from the list, which must be float
	float getFloat() {
		return _values[_getCounter++].val_float;
	}
	//! Get the next value from the list, which must be double
	double getDouble() {
		return _values[_getCounter++].val_double;
	}
	//! Get the next value from the list, which must be long double
	long double getLongDouble() {
		return _values[_getCounter++].val_longDouble;
	}
	//! number of values (possibly different types) to be communicated
	int _numValues;
	//! counter which points to the position which shall be written next
	int _setCounter;
	//! counter which points to the position which shall be read next
	int _getCounter;

	//! Array to store the values which shall be communicated
	valType* _values;

};

#endif /*COLLECTIVECOMMDUMMY_H_*/
