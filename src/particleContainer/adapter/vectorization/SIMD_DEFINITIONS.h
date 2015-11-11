/***********************************************************************************//**
 *
 * \file SIMD_DEFINITIONS.h
 *
 * \brief Contains macro definitions for the intrinsics used for the vectorization
 *
 * \author Steffen Seckler
 *
 **************************************************************************************/
#pragma once

#ifndef  SIMD_DEFINITIONS_H
#define  SIMD_DEFINITIONS_H

#ifdef IN_IDE_PARSER //just for the ide parser include the simd_types.h -- normally this is not done.
    #include "./SIMD_TYPES.h"
#endif

/*
 * Check whether the file SIMD_TYPES.hpp has been included.
 * This file (SIMD_DEFINITIONS.hpp) needs some macros to be set properly by that file (SIMD_TYPES.hpp).
 */
#ifndef  SIMD_TYPES_H
    #error "SIMD_DEFINITIONS included without SIMD_TYPES! Never include this file directly! Include it only via SIMD_TYPES!"
#endif /* defined SIMD_TYPES_H */


#if VCP_VEC_TYPE==VCP_VEC_SSE3
	static inline vcp_double_vec vcp_simd_zerov() { return _mm_setzero_pd(); }

	static inline vcp_double_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_cmplt_pd(a, b);}
	static inline vcp_double_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_cmpeq_pd(a, b);}
	static inline vcp_double_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_cmpneq_pd(a, b);}
	static inline vcp_double_vec vcp_simd_and(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_and_pd(a, b);}
	static inline vcp_double_vec vcp_simd_or(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_or_pd(a, b);}
	static inline vcp_double_vec vcp_simd_xor(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_xor_pd(a, b);}

#elif VCP_VEC_TYPE==VCP_VEC_AVX
	static inline vcp_double_vec vcp_simd_zerov() { return _mm256_setzero_pd(); }
	static inline vcp_double_vec vcp_simd_ones() { return _mm256_castsi256_pd( _mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF) ); }
	static inline vcp_double_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_cmp_pd(a, b, _CMP_LT_OS);}
	static inline vcp_double_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_cmp_pd(a, b, _CMP_EQ_OS);}
	static inline vcp_double_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_cmp_pd(a, b, _CMP_NEQ_OS);}
	static inline vcp_double_vec vcp_simd_and(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_and_pd(a, b);}
	static inline vcp_double_vec vcp_simd_or(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_or_pd(a, b);}
	static inline vcp_double_vec vcp_simd_xor(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_xor_pd(a, b);}
#endif























#endif /* SIMD_DEFINITIONS_H */
