#ifndef COLLECTIVECOMMUNICATION_H_
#define COLLECTIVECOMMUNICATION_H_

#include <mpi.h>

#include <cassert>

#include "utils/Logger.h"

/* Enable agglomerated reduce operations. This will store all values in one array and apply a
 * user defined reduce operation so that the MPI reduce operation is only called once. */
#define ENABLE_AGGLOMERATED_REDUCE 1


//! @brief This class is used to transfer several values of different types with a single command
//! @author Martin Buchholz
//!
//! At different points in the simulation process, collective communication commands are
//! necessary. One thing that always has to be done is a reduce command to get global
//! macroscopic values from the local ones. But depending on the application, other commands
//! might also be possible (also broadcast commands). The number of values to be transferred
//! from/to each proc is usually very small (about ten or even less). The main costs of
//! those communications are therefore not caused by the amount of data, but by the MPI
//! overhead caused by each call of a communication command. The purpose of this class is
//! to use a single MPI command to transfer several values of possible different types.
//! Currently supported commands are:
//! - broadcast
//! - reduce using add as reduce operation
//!
//! Currently supported datatypes are:
//! - MPI_INT
//! - MPI_UNSIGNED_LONG
//! - MPI_FLOAT
//! - MPI_DOUBLE
//! - MPI_LONG_DOUBLE
//!
//! Further datatypes and reduce operations could be added very easily.
//! A typical usage of this class could look like this:
//! @code
//!   // create object
//!   CollectiveCommunication collComm;
//!   // Set number of values to be sent
//!   collComm.init(comm, 4);
//!
//!   // store values to be sent
//!   collComm.appendInt(5);
//!   collComm.appendInt(8);
//!   collComm.appendDouble(1.3);
//!   collComm.appendUnsLong(9);
//!
//!   // execute collective communication
//!   collComm.allreduceSum();
//!
//!   // read values (IMPORTANT: same order as for storing)
//!   int i1 = collComm.getInt();
//!   int i2 = collComm.getInt();
//!   double d = collComm.getDouble();
//!   unsigned long l = collComm.getUnsLong();
//!
//!   // finalize the communication (important for deleting memory)
//!   collComm.finalize();
//! @endcode
class CollectiveCommunication {

	//! As in C++ arrays have to contain only one type of variable,
	//! but this class (see documentation of CollectiveCommunication class)
	//! shall be able to store values of different types for the transfer,
	//! this union is needed to be able to store the values of different
	//! types in one array.
	union valType {
		int val_int;
		unsigned long val_unsLong;
		float val_float;
		double val_double;
		long double val_longDouble;
	};

public:
	CollectiveCommunication() {
		_sendValues = 0;
		_recvValues = 0;
		_valuesType = MPI_DATATYPE_NULL;
	}

	virtual ~CollectiveCommunication() {
		assert(_sendValues == 0);
		assert(_recvValues == 0);
		assert( _valuesType == MPI_DATATYPE_NULL );
	}
		
	//! @brief allocate memory for the values to be sent, initialize counters
	//! @param numValues number of values that shall be communicated
	void init(MPI_Comm communicator, int numValues) {
		_communicator = communicator;
		_listOfTypes = new MPI_Datatype[numValues];
		_sendValues = new valType[numValues];
		_recvValues = new valType[numValues];
		_numValues = numValues;
		_setCounter = 0;
		_getCounter = 0;
	}

	//! @brief delete memory and MPI_Type
	void finalize() {
		_numValues = 0;
		_setCounter = 0;
		_getCounter = 0;
		delete[] _listOfTypes;
		delete[] _sendValues;
		delete[] _recvValues;
		_listOfTypes = NULL;
		_sendValues = NULL;
		_recvValues = NULL;
		assert( _valuesType == MPI_DATATYPE_NULL );
	}

	//! @brief method used by MPI to add variables of this type
	//!
	//! For collective reduce commands, an operation has to be specified
	//! which should be used to combine the values from the different processes.
	//! As with this class, special datatypes are sent, the build-in
	//! MPI operations don't work on those datatypes. Therefore, operations have
	//! to be definded which work on those datatypes. MPI allows to create own
	//! operations (type MPI_Op) by specifying a function which will be used
	//! in the reduce operation. MPI specifies the signature of such functions
	//! This methods checks from which basic datatypes the given datatype
	//! was constructed and performs an add operation for each of the basic types.
	static void add(valType *invec, valType *inoutvec, int *len, MPI_Datatype *dtype) {
		int numints;
		int numaddr;
		int numtypes;
		int combiner;

		MPI_CHECK( MPI_Type_get_envelope(*dtype, &numints, &numaddr, &numtypes, &combiner) );

		int arrayInts[numints];
		MPI_Aint arrayAddr[numaddr];
		MPI_Datatype arrayTypes[numtypes];

		MPI_CHECK( MPI_Type_get_contents(*dtype, numints, numaddr, numtypes, arrayInts, arrayAddr, arrayTypes) );

		for (int i = 0; i < numtypes; i++) {
			if (arrayTypes[i] == MPI_INT) {
				inoutvec[i].val_int += invec[i].val_int;
			}
			else if (arrayTypes[i] == MPI_UNSIGNED_LONG) {
				inoutvec[i].val_unsLong += invec[i].val_unsLong;
			}
			else if (arrayTypes[i] == MPI_FLOAT) {
				inoutvec[i].val_float += invec[i].val_float;
			}
			else if (arrayTypes[i] == MPI_DOUBLE) {
				inoutvec[i].val_double += invec[i].val_double;
			}
			else if (arrayTypes[i] == MPI_LONG_DOUBLE) {
				inoutvec[i].val_longDouble += invec[i].val_longDouble;
			}
		}
	}

	//! @brief defines a MPI datatype which can be used to transfer a CollectiveCommunication object
	//!
	//! befor this method is called, init has to be called and all values to be
	//! communicated have to be added with the append<type> methods. Then this method
	//! will construct a MPI-datatype which represents all the added values. The
	//! datatype is stored in the member variable _valuesType;
	void setMPIType() {
		int numblocks = _numValues;
		int blocklengths[numblocks];
		MPI_Aint disps[numblocks];
		int disp = 0;
		for (int i = 0; i < numblocks; i++) {
			blocklengths[i] = 1;
			disps[i] = disp;
			disp += sizeof(valType);
		}
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
	MPI_CHECK( MPI_Type_create_struct(numblocks, blocklengths, disps, _listOfTypes, &_valuesType) );
#else
	MPI_CHECK( MPI_Type_struct(numblocks, blocklengths, disps, _listOfTypes, &_valuesType) );
#endif
		MPI_CHECK( MPI_Type_commit(&_valuesType) );
	}

	//! Append an int value to the list of values to be sent
	void appendInt(int intValue) {
		_sendValues[_setCounter].val_int = intValue;
		_listOfTypes[_setCounter] = MPI_INT;
		_setCounter++;
	}

	//! Append a unsigned long value to the list of values to be sent
	void appendUnsLong(unsigned long unsLongValue) {
		_sendValues[_setCounter].val_unsLong = unsLongValue;
		_listOfTypes[_setCounter] = MPI_UNSIGNED_LONG;
		_setCounter++;
	}

	//! Append a float value to the list of values to be sent
	void appendFloat(float floatValue) {
		_sendValues[_setCounter].val_float = floatValue;
		_listOfTypes[_setCounter] = MPI_FLOAT;
		_setCounter++;
	}

	//! Append a double value to the list of values to be sent
	void appendDouble(double doubleValue) {
		_sendValues[_setCounter].val_double = doubleValue;
		_listOfTypes[_setCounter] = MPI_DOUBLE;
		_setCounter++;
	}

	//! Append a long double value to the list of values to be sent
	void appendLongDouble(double longDoubleValue) {
		_sendValues[_setCounter].val_longDouble = longDoubleValue;
		_listOfTypes[_setCounter] = MPI_LONG_DOUBLE;
		_setCounter++;
	}

	MPI_Comm getTopology() {
		return this->_communicator;
	}

	//! Get the next value from the list, which must be int
	int getInt() {
		return _recvValues[_getCounter++].val_int;
	}

	//! Get the next value from the list, which must be unsigned long
	unsigned long getUnsLong() {
		return _recvValues[_getCounter++].val_unsLong;
	}

	//! Get the next value from the list, which must be float
	float getFloat() {
		return _recvValues[_getCounter++].val_float;
	}

	//! Get the next value from the list, which must be double
	double getDouble() {
		return _recvValues[_getCounter++].val_double;
	}

	//! Get the next value from the list, which must be long double
	long double getLongDouble() {
		return _recvValues[_getCounter++].val_longDouble;
	}

	//! Broadcast all values from the process with rank 0 to all others
	void broadcast(int root = 0) {
		setMPIType();
		MPI_CHECK( MPI_Bcast(_sendValues, 1, _valuesType, root, _communicator) );
		for (int i = 0; i < _numValues; i++) {
			_recvValues[i] = _sendValues[i];
		}
		MPI_CHECK( MPI_Type_free(&_valuesType) );
	}

	//! Do Allreduce off all values with reduce operation add
	void allreduceSum() {
#if ENABLE_AGGLOMERATED_REDUCE
		setMPIType();
		MPI_Op reduceOp;
		MPI_CHECK( MPI_Op_create((MPI_User_function *) CollectiveCommunication::add, 1, &reduceOp) );
		MPI_CHECK( MPI_Allreduce(_sendValues, _recvValues, 1, _valuesType, reduceOp, _communicator) );
		MPI_CHECK( MPI_Op_free(&reduceOp) );
		MPI_CHECK( MPI_Type_free(&_valuesType) );
#else
		for( int i = 0; i < _numValues; i++ ) {
			MPI_CHECK( MPI_Allreduce( &_sendValues[i], &_recvValues[i], 1, _listOfTypes[i], MPI_SUM, _communicator ) );
		}
#endif
	}

	//! number of values (possibly different types) to be communicated
	int _numValues;
	//! counter which points to the position which shall be written next
	int _setCounter;
	//! counter which points to the position which shall be read next
	int _getCounter;

	//! Array to store the values which shall be communicated
	valType* _sendValues;
	valType* _recvValues;
	//! Array of the corresponding MPI types for the values stored in _values
	MPI_Datatype* _listOfTypes;
	//! MPI_Datatype which will be used in the communication command and which
	//! represents all values
	MPI_Datatype _valuesType;

	//! Communicater to be used by the communication commands
	MPI_Comm _communicator;

};

#endif /* COLLECTIVECOMMUNICATION_H_ */
