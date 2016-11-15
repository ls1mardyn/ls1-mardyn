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

#define vcp_inline inline __attribute__((always_inline))

#ifdef IN_IDE_PARSER //just for the ide parser include the simd_types.h -- normally this is not done.
    #include "./SIMD_TYPES.h"
#endif

#include "DoubleVec.h"
#include "math.h"

using namespace vcp;

/*
 * Check whether the file SIMD_TYPES.hpp has been included.
 * This file (SIMD_DEFINITIONS.hpp) needs some macros to be set properly by that file (SIMD_TYPES.hpp).
 */
#ifndef  SIMD_TYPES_H
    #error "SIMD_DEFINITIONS included without SIMD_TYPES! Never include this file directly! Include it only via SIMD_TYPES!"
#endif /* defined SIMD_TYPES_H */

#if VCP_VEC_TYPE==VCP_NOVEC
	static vcp_inline vcp_mask_vec vcp_simd_getInitMask(const size_t& /*i*/){
		return true;
	}
	static vcp_inline vcp_mask_vec vcp_simd_getRemainderMask(const size_t& /*size*/){
		return false;
	}

#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	static vcp_inline vcp_mask_vec vcp_simd_getInitMask(const size_t& i){
		switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
			case 0: return _mm_set_epi32(~0, ~0, ~0, ~0);
			default: return _mm_set_epi32(~0, ~0, 0, 0);
		}
	}
	static vcp_inline vcp_mask_vec vcp_simd_getRemainderMask(const size_t& size) {
		switch (size & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
			case 0: return _mm_set_epi32(0, 0, 0, 0);
			default: return _mm_set_epi32(0, 0, ~0, ~0);
		}
	}
#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2
	static vcp_inline vcp_mask_vec vcp_simd_getInitMask(const size_t& i){
		switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
			case 0: return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0);
			case 1: return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, 0, 0);
			case 2: return _mm256_set_epi32(~0, ~0, ~0, ~0, 0, 0, 0, 0);
			default: return _mm256_set_epi32(~0, ~0, 0, 0, 0, 0, 0, 0);
		}
	}
	static vcp_inline vcp_mask_vec vcp_simd_getRemainderMask(const size_t& size) {
		switch (size & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
			case 0: return MaskVec::zero();
			case 1: return _mm256_set_epi32(0, 0, 0, 0, 0, 0, ~0, ~0);
			case 2: return _mm256_set_epi32(0, 0, 0, 0, ~0, ~0, ~0, ~0);
			default: return _mm256_set_epi32(0, 0, ~0, ~0, ~0, ~0, ~0, ~0);
		}
	}
#elif VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or\
	  VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER

	static vcp_inline vcp_mask_vec vcp_simd_getInitMask(const size_t& i){
		static const vcp_mask_vec possibleInitJMasks[VCP_VEC_SIZE] = { 0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE0, 0xC0, 0x80 };
		return possibleInitJMasks[i & static_cast<size_t>(VCP_VEC_SIZE_M1)];
	}

	static vcp_inline vcp_mask_vec vcp_simd_getRemainderMask(const size_t& size) {
		static const vcp_mask_vec possibleRemainderJMasks[VCP_VEC_SIZE] = { 0x00, 0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3F, 0x7F };
		return possibleRemainderJMasks[size & static_cast<size_t>(VCP_VEC_SIZE_M1)];
	}
#endif

/**
 * if num is not divisible by VCP_VEC_SIZE finds the next bigger multiple of VCP_VEC_SIZE, otherwise returns num.
 * @param num
 * @return
 */
template<class T>
static vcp_inline T vcp_ceil_to_vec_size(const T& num){
	return (num + static_cast<T>(VCP_VEC_SIZE_M1)) & (~static_cast<T>(VCP_VEC_SIZE_M1));
}

/**
 * if num is not divisible by VCP_VEC_SIZE finds the next smaller multiple of VCP_VEC_SIZE, otherwise returns num.
 * @param num
 * @return
 */
template<class T>
static vcp_inline T vcp_floor_to_vec_size(const T& num){
	return num & (~static_cast<T>(VCP_VEC_SIZE_M1));
}

#endif /* SIMD_DEFINITIONS_H */
