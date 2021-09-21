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

#include <type_traits>
#include <cstdint>

// definitions for precision:
#if defined(VCP_SPSP) or \
	defined(VCP_SPDP) or \
	defined(VCP_DPDP)
	#error conflicting precision definitions
#endif

#define VCP_SPSP 0
#define VCP_SPDP 1
#define VCP_DPDP 2

#if defined(MARDYN_SPSP)
	#define VCP_PREC VCP_SPSP
#elif defined(MARDYN_SPDP)
	#define VCP_PREC VCP_SPDP
#elif defined(MARDYN_DPDP)
	#define VCP_PREC VCP_DPDP
#else
	#error Precision macro not defined.
#endif

// The following error should NEVER occur, since it signalizes, that the macros, used by THIS translation unit are defined anywhere else in the program.
#if defined(VCP_VEC_TYPE) || defined(VCP_NOVEC) || defined(VCP_VEC_SSE3) || defined(VCP_VEC_AVX) ||         \
	defined(VCP_VEC_AVX2) || defined(VCP_VEC_KNL) || defined(VCP_VEC_KNL_GATHER) || defined(VCP_VEC_SVE) || \
	defined(VCP_VEC_SVE_GATHER)
#error conflicting macro definitions
#endif


#define VCP_NOVEC 0
#define VCP_VEC_SSE3 1
#define VCP_VEC_AVX 2
#define VCP_VEC_AVX2 3
#define VCP_VEC_KNL 6
#define VCP_VEC_KNL_GATHER 7
#define VCP_VEC_AVX512F 8
#define VCP_VEC_AVX512F_GATHER 9
#define VCP_VEC_SVE 10
#define VCP_VEC_SVE_GATHER 11

#define VCP_VEC_W__64 0
#define VCP_VEC_W_128 1
#define VCP_VEC_W_256 2
#define VCP_VEC_W_512 3

#if defined(__AVX2__) && not defined(__FMA__)  // fma should always be existent alongside avx2!!!
	#error AVX2 enabled, but no FMA found. Please enable fma to use avx2.
#endif

// define symbols for vectorization
#if defined(__AVX512F__) && defined(__AVX512ER__)
	#if defined(__VCP_GATHER__)
		#define VCP_VEC_TYPE VCP_VEC_KNL_GATHER
	#else
		#define VCP_VEC_TYPE VCP_VEC_KNL
	#endif
#elif defined(__AVX512F__)
	#if defined(__VCP_GATHER__)
		#define VCP_VEC_TYPE VCP_VEC_AVX512F_GATHER
	#else
		#define VCP_VEC_TYPE VCP_VEC_AVX512F
	#endif
#elif defined(__AVX2__) && defined(__FMA__)
	#define VCP_VEC_TYPE VCP_VEC_AVX2
#elif defined(__AVX__) && not defined(AVX128)
	#define VCP_VEC_TYPE VCP_VEC_AVX
#elif defined(__AVX__) && defined(AVX128)
	#define VCP_VEC_TYPE VCP_VEC_SSE3
#elif defined(__SSE3__)
	#define VCP_VEC_TYPE VCP_VEC_SSE3
#elif defined(__ARM_FEATURE_SVE)
	#if defined(__VCP_GATHER__)
		#define VCP_VEC_TYPE VCP_VEC_SVE_GATHER
	#else
		#define VCP_VEC_TYPE VCP_VEC_SVE
	#endif
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
#if VCP_VEC_TYPE==VCP_NOVEC
	// no file to include
#elif VCP_VEC_TYPE==VCP_VEC_SVE
	#include "arm_sve.h"
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	#include "pmmintrin.h"
#else
	// all others need immintrin.h
	#include "immintrin.h"
#endif

// define necessary types

typedef std::conditional<VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP, float, double>::type vcp_real_calc;
typedef std::conditional<VCP_PREC == VCP_SPSP                        , float, double>::type vcp_real_accum;
typedef std::conditional<VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP, uint32_t, uint64_t>::type vcp_ljc_id_t;

typedef int countertype32;//int is 4Byte almost everywhere... replace with __int32 if problems occur

#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.

	#define VCP_VEC_WIDTH VCP_VEC_W__64

	typedef uint8_t vcp_mask_vec;
	typedef uint8_t vcp_mask_single;
	typedef vcp_mask_vec vcp_lookupOrMask_vec;
	typedef vcp_mask_single vcp_lookupOrMask_single;

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3

	#define VCP_VEC_WIDTH VCP_VEC_W_128

	typedef std::conditional<VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP, uint32_t, uint64_t>::type vcp_mask_single;

	typedef __m128i vcp_mask_vec;
	typedef vcp_mask_vec vcp_lookupOrMask_vec;
	typedef vcp_mask_single vcp_lookupOrMask_single;

#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2//avx, avx2

	#define VCP_VEC_WIDTH VCP_VEC_W_256

	typedef std::conditional<VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP, uint32_t, uint64_t>::type vcp_mask_single;

	typedef __m256i vcp_mask_vec;
	typedef vcp_mask_vec vcp_lookupOrMask_vec;
	typedef vcp_mask_single vcp_lookupOrMask_single;

#elif VCP_VEC_TYPE==VCP_VEC_KNL or \
	  VCP_VEC_TYPE==VCP_VEC_KNL_GATHER or \
	  VCP_VEC_TYPE==VCP_VEC_AVX512F or \
	  VCP_VEC_TYPE==VCP_VEC_AVX512F_GATHER

	#define VCP_VEC_WIDTH VCP_VEC_W_512

	// these can't be put in std::conditional, because it would just be too nice. warnings.
	#if VCP_PREC==VCP_SPSP or VCP_PREC==VCP_SPDP
		typedef __mmask16 vcp_mask_vec;
		typedef __mmask16 vcp_mask_single;
	#else // VCP_PREC==VCP_DPDP
		typedef __mmask8 vcp_mask_vec;
		typedef __mmask8 vcp_mask_single;
	#endif

	#if VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_AVX512F
		typedef vcp_mask_vec vcp_lookupOrMask_vec;
		typedef vcp_mask_single vcp_lookupOrMask_single;

	#else  // VCP_VEC_TYPE==VCP_VEC_KNL_GATHER or VCP_VEC_AVX512F_GATHER
		typedef __m512i vcp_lookupOrMask_vec;
		typedef countertype32 vcp_lookupOrMask_single;
	#endif
#elif VCP_VEC_TYPE==VCP_VEC_SVE
	#ifndef SVE_VEC_WIDTH
		#error VCP_VEC_TYPE is VCP_VEC_SVE, but SVE_VEC_WIDTH is not set!
	#endif
	// TODO: there are lots of problems with this!
	#define VCP_VEC_WIDTH VCP_VEC_SVE

	// TODO! vcp_mask_single is probably wrong, could be bool?
	typedef std::conditional<VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP, uint32_t, uint64_t>::type vcp_mask_single;

	typedef svbool_t vcp_mask_vec;
	typedef vcp_mask_vec vcp_lookupOrMask_vec;
	typedef vcp_mask_single vcp_lookupOrMask_single;
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
