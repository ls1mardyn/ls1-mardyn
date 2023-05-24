/*
 * DynAlloc.h
 *
 *  Created on: 09.11.2016
 *      Author: mheinen
 */

#ifndef DYNALLOC_H_
#define DYNALLOC_H_

#include <iostream>
#include <cstdint>

// leak save dynamic memory allocation
inline void AllocateUnsLongArray(unsigned long* &ptr, const unsigned int& nSize)
{
	if(NULL != ptr)
		delete[] ptr;
	ptr = new unsigned long[nSize];
}

inline void AllocateDoubleArray(double* &ptr, const unsigned int& nSize)
{
	if(NULL != ptr)
		delete[] ptr;
	ptr = new double[nSize];
}

inline void AllocateUint8Array(uint8_t* &ptr, const uint32_t& nSize)
{
	if(NULL != ptr)
		delete[] ptr;
	ptr = new uint8_t[nSize];
}

inline void AllocateUint16Array(uint16_t* &ptr, const uint32_t& nSize)
{
	if(NULL != ptr)
		delete[] ptr;
	ptr = new uint16_t[nSize];
}

inline void AllocateUint32Array(uint32_t* &ptr, const uint32_t& nSize)
{
	if(NULL != ptr)
		delete[] ptr;
	ptr = new uint32_t[nSize];
}

inline void AllocateUint64Array(uint64_t* &ptr, const uint32_t& nSize)
{
	if(NULL != ptr)
		delete[] ptr;
	ptr = new uint64_t[nSize];
}


#endif /* DYNALLOC_H_ */
