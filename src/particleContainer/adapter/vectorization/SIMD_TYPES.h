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
#define VCP_VEC_KNC 4
#define VCP_VEC_KNC_GATHER 5
#define VCP_VEC_KNL 6
#define VCP_VEC_KNL_GATHER 7

#define VCP_VEC_W__64 0
#define VCP_VEC_W_128 1
#define VCP_VEC_W_256 2
#define VCP_VEC_W_512 3


typedef int countertype32;//int is 4Byte almost everywhere... replace with __int32 if problems occur


#if defined(__AVX2__) && not defined(__FMA__)//fma should always be existent alongside avx2!!!
	#warn AVX2 enabled, but no FMA found. Please enable fma to use avx2.
#endif

// define symbols for vectorization
#if defined(__AVX512F__)
	#if not defined(__AVX512ER__)
	#error currently only KNL is supported, not SKX
	#endif
	#if defined(__VCP_GATHER__)
		#define VCP_VEC_TYPE VCP_VEC_KNL_GATHER
	#else
		#define VCP_VEC_TYPE VCP_VEC_KNL
	#endif
#elif defined(__MIC__)
	#if defined(__VCP_GATHER__)
		#define VCP_VEC_TYPE VCP_VEC_KNC_GATHER
	#else
		#define VCP_VEC_TYPE VCP_VEC_KNC
	#endif
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
#if VCP_VEC_TYPE==VCP_VEC_AVX or \
	VCP_VEC_TYPE==VCP_VEC_AVX2 or \
	VCP_VEC_TYPE==VCP_VEC_KNC or \
	VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or\
	VCP_VEC_TYPE==VCP_VEC_KNL or \
	VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
	#include "immintrin.h"
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	#include "pmmintrin.h"
#endif



// define necessary types

#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	typedef double vcp_double_vec;
	#define VCP_VEC_SIZE 1u
	#define VCP_VEC_SIZE_M1 0u
	#define VCP_VEC_WIDTH VCP_VEC_W__64

	typedef bool vcp_mask_vec;
	typedef bool vcp_mask_single;
	#define VCP_INDICES_PER_LOOKUP_SINGLE 1
	#define VCP_INDICES_PER_LOOKUP_SINGLE_M1 0
	typedef vcp_mask_vec vcp_lookupOrMask_vec;
	typedef vcp_mask_single vcp_lookupOrMask_single;

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	typedef __m128d vcp_double_vec;
	#define VCP_VEC_SIZE 2u
	#define VCP_VEC_SIZE_M1 1u
	#define VCP_VEC_WIDTH VCP_VEC_W_128

	typedef __m128i vcp_mask_vec;
	typedef unsigned long vcp_mask_single;
	#define VCP_INDICES_PER_LOOKUP_SINGLE 1
	#define VCP_INDICES_PER_LOOKUP_SINGLE_M1 0
	typedef vcp_mask_vec vcp_lookupOrMask_vec;
	typedef vcp_mask_single vcp_lookupOrMask_single;

#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2//avx, avx2
	typedef __m256d vcp_double_vec;
	#define VCP_VEC_SIZE 4u
	#define VCP_VEC_SIZE_M1 3u
	#define VCP_VEC_WIDTH VCP_VEC_W_256

	typedef __m256i vcp_mask_vec;
	typedef unsigned long vcp_mask_single;
	#define VCP_INDICES_PER_LOOKUP_SINGLE 1
	#define VCP_INDICES_PER_LOOKUP_SINGLE_M1 0
	typedef vcp_mask_vec vcp_lookupOrMask_vec;
	typedef vcp_mask_single vcp_lookupOrMask_single;

#elif VCP_VEC_TYPE==VCP_VEC_KNC or \
	  VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or \
	  VCP_VEC_TYPE==VCP_VEC_KNL or \
	  VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
	typedef __m512d vcp_double_vec;
	#define VCP_VEC_SIZE 8u
	#define VCP_VEC_SIZE_M1 7u
	#define VCP_VEC_WIDTH VCP_VEC_W_512

	typedef __mmask8 vcp_mask_vec;
	typedef __mmask8 vcp_mask_single;
	#if VCP_VEC_TYPE==VCP_VEC_KNC or\
		VCP_VEC_TYPE==VCP_VEC_KNL
		#define VCP_INDICES_PER_LOOKUP_SINGLE 8
		#define VCP_INDICES_PER_LOOKUP_SINGLE_M1 7
		typedef __mmask8 vcp_lookupOrMask_vec;
		typedef __mmask8 vcp_lookupOrMask_single;
	#else  // VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
		#define VCP_INDICES_PER_LOOKUP_SINGLE 1
		#define VCP_INDICES_PER_LOOKUP_SINGLE_M1 0
		typedef __m512i vcp_lookupOrMask_vec;
		typedef countertype32 vcp_lookupOrMask_single;
	#endif

#endif







/*
 * Control macros are set
 *
 * Include the function macros
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "SIMD_DEFINITIONS.h"

#endif /* SIMD_TYPES_H */
