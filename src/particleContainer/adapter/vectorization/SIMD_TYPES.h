/***********************************************************************************//**
 *
 * \file SIMD_TYPES.h
 *
 * \brief Defines the length of the vectors and the corresponding functions
 *
 * \author Steffen Seckler
 *
 **************************************************************************************/
#pragma once
#ifndef  SIMD_TYPES_H
#define  SIMD_TYPES_H


// The following error should NEVER occur, since it signalizes, that the macros, used by THIS translation unit are defined anywhere else in the program.
#if defined(VCP_VEC_TYPE) || defined(VCP_NOVEC) || defined(VCP_VEC_SSE3) || defined(VCP_VEC_AVX) || defined(VCP_VEC_AVX2)
	#error conflicting macro definitions
#endif


#define VCP_NOVEC 0
#define VCP_VEC_SSE3 1
#define VCP_VEC_AVX 2
#define VCP_VEC_AVX2 3
#define VCP_VEC_MIC 4

//AVX2: __AVX2__ (icc,gcc)
//fma: __FMA__ (icc,gcc) <- always existent for avx2 architectures, but sometimes has to be enabled.


#if defined(__AVX2__) && not defined(__FMA__)//fma should always be existent alongside avx2!!!
	#warn AVX2 enabled, but no FMA found. Please enable fma to use avx2.
#endif

// define symbols for vectorization
#if defined(__MIC__)
	#define VCP_VEC_TYPE VCP_VEC_MIC
#elif defined(__AVX2__) && defined(__FMA__)
	#define VCP_VEC_TYPE VCP_VEC_AVX2
#elif defined(__AVX__) && not defined(AVX128)
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
#if VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2 or VCP_VEC_TYPE==VCP_VEC_MIC
	#include "immintrin.h"
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	#include "pmmintrin.h"
#endif



// define necessary types

#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	typedef double vcp_double_vec;
	typedef unsigned long vcp_doublesizedmask_vec;
	#define VCP_VEC_SIZE 1u
	#define VCP_VEC_SIZE_M1 0u

	typedef unsigned long vcp_mask_vec;

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	typedef __m128d vcp_double_vec;
	typedef vcp_double_vec vcp_doublesizedmask_vec;
	#define VCP_VEC_SIZE 2u
	#define VCP_VEC_SIZE_M1 1u

	typedef __m128i vcp_mask_vec;

#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2//avx, avx2
	typedef __m256d vcp_double_vec;
	typedef vcp_double_vec vcp_doublesizedmask_vec;
	#define VCP_VEC_SIZE 4u
	#define VCP_VEC_SIZE_M1 3u

	typedef __m256i vcp_mask_vec;

#elif VCP_VEC_TYPE==VCP_VEC_MIC//mic
	typedef __m512d vcp_double_vec;
	typedef vcp_double_vec vcp_doublesizedmask_vec;
	#define VCP_VEC_SIZE 8u
	#define VCP_VEC_SIZE_M1 7u

	typedef __mmask8 vcp_mask_vec;

#endif












/*
 * Control macros are set
 *
 * Include the function macros
 */
#include "SIMD_DEFINITIONS.h"


#endif /* SIMD_TYPES_H */
