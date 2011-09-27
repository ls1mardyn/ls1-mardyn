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

#ifdef CONFIG_CUDA_FLOAT_UNSORTED
#	define CONFIG_NAME "cuda_float_unsorted"
#endif

#ifdef CONFIG_CUDA_FLOAT_SORTED_HWCACHEONLY
#       define CONFIG_NAME "cuda_float_sorted"

#       define CUDA_HW_CACHE_ONLY

#       define CUDA_SORT_CELLS_BY_COMPONENTTYPE

#endif

#ifdef CONFIG_CUDA_FLOAT_SORTED_WBDP
#       define CONFIG_NAME "cuda_float_sorted_wbdp"

#       define CUDA_SORT_CELLS_BY_COMPONENTTYPE

#       define CUDA_HW_CACHE_ONLY
#       define CUDA_WARP_BLOCK_CELL_PROCESSOR
#endif

#ifdef CONFIG_CUDA_DOUBLE_UNSORTED
#	define CONFIG_NAME "cuda_double_unsorted"
#	define CUDA_DOUBLE_MODE
#endif

#ifdef CONFIG_CUDA_DOUBLE_UNSORTED_WBDP
#	define CONFIG_NAME "cuda_double_unsorted_wbdp"
#	define CUDA_DOUBLE_MODE

#	define CUDA_HW_CACHE_ONLY
#	define CUDA_WARP_BLOCK_CELL_PROCESSOR
#endif

#ifdef CONFIG_CUDA_DOUBLE_SORTED_WBDP
#       define CONFIG_NAME "cuda_double_sorted_wbdp"
#       define CUDA_DOUBLE_MODE

#       define CUDA_SORT_CELLS_BY_COMPONENTTYPE

#       define CUDA_HW_CACHE_ONLY
#       define CUDA_WARP_BLOCK_CELL_PROCESSOR
#endif

#ifdef CONFIG_CUDA_DOUBLE_SORTED
#	define CONFIG_NAME "cuda_double_sorted"
#	define CUDA_DOUBLE_MODE
#	define CUDA_SORT_CELLS_BY_COMPONENTTYPE
#endif

#ifdef CONFIG_CUDA_DOUBLE_UNSORTED_HWCACHEONLY
#       define CONFIG_NAME "cuda_double_unsorted_hwcacheonly"

#       define CUDA_DOUBLE_MODE

#       define CUDA_HW_CACHE_ONLY
#endif

#ifdef CONFIG_CUDA_DOUBLE_SORTED_HWCACHEONLY
#	define CONFIG_NAME "cuda_double_sorted_hwcacheonly"

#	define CUDA_DOUBLE_MODE
#	define CUDA_SORT_CELLS_BY_COMPONENTTYPE

#	define CUDA_HW_CACHE_ONLY
#endif

#ifdef CONFIG_CUDA_DOUBLE_UNSORTED_UNPACKED_STORAGE
#       define CONFIG_NAME "cuda_double_unsorted_unpacked_storage"

#       define CUDA_DOUBLE_MODE

#       define CUDA_UNPACKED_STORAGE
#endif

#ifdef CONFIG_CUDA_DOUBLE_UNSORTED_WBDP_UNPACKED_STORAGE
#       define CONFIG_NAME "cuda_double_unsorted_wbdp_unpacked_storage"

#       define CUDA_DOUBLE_MODE

#       define CUDA_HW_CACHE_ONLY

#       define CUDA_UNPACKED_STORAGE
#       define CUDA_WARP_BLOCK_CELL_PROCESSOR
#endif

#ifdef CONFIG_CUDA_DOUBLE_SORTED_WBDP
#       define CONFIG_NAME "cuda_double_sorted_wbdp"

#       define CUDA_DOUBLE_MODE

#       define CUDA_HW_CACHE_ONLY

#       define CUDA_SORT_CELLS_BY_COMPONENTTYPE

#       define CUDA_WARP_BLOCK_CELL_PROCESSOR
#endif

#ifdef CONFIG_CUDA_DOUBLE_SORTED_WBDP_WITH_CACHE
#       define CONFIG_NAME "cuda_double_sorted_wbdp_with_cache"

#       define CUDA_DOUBLE_MODE

#       define CUDA_SORT_CELLS_BY_COMPONENTTYPE

#       define CUDA_WARP_BLOCK_CELL_PROCESSOR
#endif

#ifndef CONFIG_NAME
#error no CONFIG_* specified!
#endif

//#define CUDA_WARP_BLOCK_CELL_PROCESSOR

//#define CUDA_UNPACKED_STORAGE

//#define CUDA_HW_CACHE_ONLY

//#define REFERENCE_IMPLEMENTATION
//#define TEST_QUATERNION_MATRIX_CONVERSION
#define COMPARE_TO_CPU
//#define USE_BEHAVIOR_PROBE
//#define DUMP_CELL_LENGTHS
//#define DEBUG_COMPONENT_DESCRIPTORS

//#define MAX_NUM_WARPS 6

//#define MAX_NUM_COMPONENTS 2

//#define MAX_NUM_LJCENTERS 3
//#define MAX_NUM_DIPOLES 3
//#define MAX_NUM_CHARGES 0

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

//#define __restrict__

// TODO: move this include into the referencing header files
#include "sharedDecls.h"

#endif /* CONFIG_H_ */
