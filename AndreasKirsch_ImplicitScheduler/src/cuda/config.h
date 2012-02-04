/*
 * config.h
 *
 *  Created on: Jun 14, 2011
 *      Author: andreas
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include "cutil_double_math.h"

#define MAX_BLOCK_SIZE 1024
#define QUATERNION_BLOCK_SIZE 512

// created by the makefile to avoid having to call make clean all the time
#include "make_config.h"

#define BENCHMARKING

//#define CONFIG_CUDA_DOUBLE_SORTED

#ifdef CONFIG_NO_CUDA
#	define CONFIG_NAME "no_cuda"
#	define NO_CUDA
#endif

#ifdef CONFIG_CUDA_DOUBLE_REFERENCE
#	define CONFIG_NAME "cuda_double_ref"

#	define CUDA_DOUBLE_MODE
#	define REFERENCE_IMPLEMENTATION
#endif

#ifdef CONFIG_CUDA_FLOAT
#	define CONFIG_NAME "cuda_float"
#endif

#ifdef CONFIG_CUDA_FLOAT_WBCP
#	define CONFIG_NAME "cuda_float_wbcp"

#	define CUDA_WARP_BLOCK_CELL_PROCESSOR
#endif

#ifdef CONFIG_CUDA_DOUBLE_WBCP
#	define CONFIG_NAME "cuda_double_wbcp"
#	define CUDA_DOUBLE_MODE

#	define CUDA_WARP_BLOCK_CELL_PROCESSOR
#endif

#ifdef CONFIG_CUDA_DOUBLE
#	define CONFIG_NAME "cuda_double"
#	define CUDA_DOUBLE_MODE
#endif

#ifndef CONFIG_NAME
#error no CONFIG_* specified!
#endif

// uncomment this to measure the GPU error compared to CPU results
#define COMPARE_TO_CPU

// CPU/GPU quaternion calculation comparison
//#define TEST_QUATERNION_MATRIX_CONVERSION

//#define USE_BEHAVIOR_PROBE

// dump all cell sizes (for creating histograms)
//#define DUMP_CELL_LENGTHS

// dump the component descriptors
//#define DEBUG_COMPONENT_DESCRIPTORS

//#define MAX_NUM_WARPS 6

//#define MAX_NUM_COMPONENTS 2

//#define MAX_NUM_LJCENTERS 3
//#define MAX_NUM_DIPOLES 3
//#define MAX_NUM_CHARGES 0
//#define MAX_NUM_QUADRUPOLES 0

#ifndef REFERENCE_IMPLEMENTATION
#	define WARP_SIZE 32u
#	define NUM_WARPS MAX_NUM_WARPS
#else
#	define WARP_SIZE 1u
#	define NUM_WARPS 1
#endif

#ifndef MAX_REGISTER_COUNT
#	define MAX_REGISTER_COUNT 63
#endif

#define BLOCK_SIZE uint(WARP_SIZE*NUM_WARPS)

#ifdef CUDA_DOUBLE_MODE
	typedef double floatType;
	typedef double3 floatType3;

#	define make_floatType3 make_double3
#else
	typedef float floatType;
	typedef float3 floatType3;

#	define make_floatType3 make_float3
#endif

// TODO: move this include into the referencing header files
#include "sharedDecls.h"

#endif /* CONFIG_H_ */
