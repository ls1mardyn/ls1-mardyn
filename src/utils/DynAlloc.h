/*
 * DynAlloc.h
 *
 *  Created on: 09.11.2016
 *      Author: mheinen
 */

#ifndef DYNALLOC_H_
#define DYNALLOC_H_

#include <iostream>
using namespace std;

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


#endif /* DYNALLOC_H_ */
