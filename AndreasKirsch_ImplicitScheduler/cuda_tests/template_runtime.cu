/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

/* Template project which demonstrates the basics on how to setup a project 
* example application, doesn't use cutil library.
*/

#include <stdio.h>
#include <string.h>
#include <iostream>

using namespace std;

bool g_bQATest = false;

#ifdef _WIN32
   #define STRCASECMP  _stricmp
   #define STRNCASECMP _strnicmp
#else
   #define STRCASECMP  strcasecmp
   #define STRNCASECMP strncasecmp
#endif

#define ASSERT(x, msg, retcode) \
    if (!(x)) \
    { \
        cout << msg << " " << __FILE__ << ":" << __LINE__ << endl; \
        return retcode; \
    }

inline cudaError cutilDeviceSynchronize()
{
	return cudaDeviceSynchronize();
}

inline void cutilDeviceReset()
{
	cudaDeviceReset();
}

extern void testGlobalSemaphore();
extern void testSingleLock();

int main(int argc, char **argv)
{
    testGlobalSemaphore();
    cutilDeviceReset();

    testSingleLock();
    cutilDeviceReset();
}
