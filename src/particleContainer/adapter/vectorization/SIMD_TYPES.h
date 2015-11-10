/***********************************************************************************//**
 *
 * \file SIMD_TYPES.h
 *
 * \brief Defines the length of the vectors and the corresponding functions
 *
 * \author Steffen Seckler (seckler), seckler AT in.tum.de
 *
 **************************************************************************************/
#pragma once
#ifndef  SIMD_TYPES_H
#define  SIMD_TYPES_H


// The following error should NEVER occur, since it signalizes, that the macros, used by THIS translation unit are defined anywhere else in the program.
#if defined(VCP_VEC_TYPE) || defined(VCP_NOVEC) || defined(VCP_VEC_SSE3) || defined(VCP_VEC_AVX)
	#error conflicting macro definitions
#endif


#define VCP_NOVEC 0
#define VCP_VEC_SSE3 1
#define VCP_VEC_AVX 2
// define symbols for vectorization
#if defined(__AVX__) && not defined(AVX128)
	#define VCP_VEC_TYPE VCP_VEC_AVX
#elif defined(__AVX__) && defined(AVX128)
	#define VCP_VEC_TYPE VCP_VEC_SSE3
#elif defined(__SSE3__)
	#define VCP_VEC_TYPE VCP_VEC_SSE3
#else
	#define VCP_VEC_TYPE VCP_NOVEC
#endif

#ifdef NOVEC
	#ifdef VCP_VEC_TYPE
		#warn Multiple vectorization methods specified. Will not use vectorization at all!
		#undef VCP_VEC_TYPE
	#endif
	#define VCP_VEC_TYPE VCP_NOVEC
#endif

// Include necessary files if we vectorize.
#if VCP_VEC_TYPE==VCP_VEC_AVX
	#include "immintrin.h"
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	#include "pmmintrin.h"
#endif




















































































/*
 * Control macros are set
 *
 * Include the function macros
 */
#include "SIMD_DEFINITIONS.h"
#endif /* SIMD_TYPES_H */
