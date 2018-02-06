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

#if defined(__GNUC__)
#define vcp_inline inline __attribute__((always_inline))
#else
#define vcp_inline inline
#endif

#ifdef IN_IDE_PARSER //just for the ide parser include the simd_types.h -- normally this is not done.
    #include "./SIMD_TYPES.h"
#endif

#include "RealVec.h"

namespace vcp {
	typedef MaskVec<vcp_real_calc> MaskCalcVec;
	typedef RealVec<vcp_real_calc> RealCalcVec;
}

// use constexpr instead of conditional compilation to death:

constexpr size_t VCP_VEC_SIZE = sizeof(vcp::RealCalcVec) / sizeof(vcp_real_calc);
constexpr size_t VCP_VEC_SIZE_M1 = VCP_VEC_SIZE - 1u;

constexpr size_t VCP_INDICES_PER_LOOKUP_SINGLE = (VCP_VEC_TYPE != VCP_VEC_KNL) ? 1u : VCP_VEC_SIZE;
constexpr size_t VCP_INDICES_PER_LOOKUP_SINGLE_M1 = (VCP_VEC_TYPE != VCP_VEC_KNL) ? 0u : VCP_VEC_SIZE_M1;

constexpr size_t VCP_ALIGNMENT = (VCP_VEC_TYPE != VCP_NOVEC) ? sizeof(vcp::RealCalcVec) : 8u;

#include <cmath>
#include "sys/types.h"


using namespace vcp;

/*
 * Check whether the file SIMD_TYPES.hpp has been included.
 * This file (SIMD_DEFINITIONS.hpp) needs some macros to be set properly by that file (SIMD_TYPES.hpp).
 */
#ifndef  SIMD_TYPES_H
    #error "SIMD_DEFINITIONS included without SIMD_TYPES! Never include this file directly! Include it only via SIMD_TYPES!"
#endif /* defined SIMD_TYPES_H */

#if VCP_VEC_TYPE==VCP_NOVEC
	static vcp_inline MaskCalcVec vcp_simd_getInitMask(const size_t& /*i*/){
		return true;
	}
	static vcp_inline MaskCalcVec vcp_simd_getRemainderMask(const size_t& /*size*/){
		return false;
	}

#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		static vcp_inline MaskCalcVec vcp_simd_getInitMask(const size_t& i){
			switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0:  return _mm_set_epi32(~0, ~0, ~0, ~0);
				case 1:  return _mm_set_epi32(~0, ~0, ~0,  0);
				case 2:  return _mm_set_epi32(~0, ~0,  0,  0);
				default: return _mm_set_epi32(~0,  0,  0,  0);
			}
		}
		static vcp_inline MaskCalcVec vcp_simd_getRemainderMask(const size_t& size) {
			switch (size & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0:  return _mm_set_epi32(0,  0,  0,  0);
				case 1:  return _mm_set_epi32(0,  0,  0, ~0);
				case 2:  return _mm_set_epi32(0,  0, ~0, ~0);
				default: return _mm_set_epi32(0, ~0, ~0, ~0);
			}
		}
	#else /*VCP_DPDP*/
		static vcp_inline MaskCalcVec vcp_simd_getInitMask(const size_t& i){
			switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0: return _mm_set_epi32(~0, ~0, ~0, ~0);
				default: return _mm_set_epi32(~0, ~0, 0, 0);
			}
		}
		static vcp_inline MaskCalcVec vcp_simd_getRemainderMask(const size_t& size) {
			switch (size & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0: return _mm_set_epi32(0, 0, 0, 0);
				default: return _mm_set_epi32(0, 0, ~0, ~0);
			}
		}
	#endif
#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2
	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		static vcp_inline MaskCalcVec vcp_simd_getInitMask(const size_t& i){
			switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0:  return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0);
				case 1:  return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0,  0);
				case 2:  return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0,  0,  0);
				case 3:  return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0,  0,  0,  0);
				case 4:  return _mm256_set_epi32(~0, ~0, ~0, ~0,  0,  0,  0,  0);
				case 5:  return _mm256_set_epi32(~0, ~0, ~0,  0,  0,  0,  0,  0);
				case 6:  return _mm256_set_epi32(~0, ~0,  0,  0,  0,  0,  0,  0);
				default: return _mm256_set_epi32(~0,  0,  0,  0,  0,  0,  0,  0);
			}
		}
		static vcp_inline MaskCalcVec vcp_simd_getRemainderMask(const size_t& size) {
			switch (size & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0:  return _mm256_set_epi32( 0,  0,  0,  0,  0,  0,  0,  0);
				case 1:  return _mm256_set_epi32( 0,  0,  0,  0,  0,  0,  0, ~0);
				case 2:  return _mm256_set_epi32( 0,  0,  0,  0,  0,  0, ~0, ~0);
				case 3:  return _mm256_set_epi32( 0,  0,  0,  0,  0, ~0, ~0, ~0);
				case 4:  return _mm256_set_epi32( 0,  0,  0,  0, ~0, ~0, ~0, ~0);
				case 5:  return _mm256_set_epi32( 0,  0,  0, ~0, ~0, ~0, ~0, ~0);
				case 6:  return _mm256_set_epi32( 0,  0, ~0, ~0, ~0, ~0, ~0, ~0);
				default: return _mm256_set_epi32( 0, ~0, ~0, ~0, ~0, ~0, ~0, ~0);
			}
		}
	#else /* VCP_DPDP */
		static vcp_inline MaskCalcVec vcp_simd_getInitMask(const size_t& i){
			switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0: return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0);
				case 1: return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, 0, 0);
				case 2: return _mm256_set_epi32(~0, ~0, ~0, ~0, 0, 0, 0, 0);
				default: return _mm256_set_epi32(~0, ~0, 0, 0, 0, 0, 0, 0);
			}
		}
		static vcp_inline MaskCalcVec vcp_simd_getRemainderMask(const size_t& size) {
			switch (size & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0: return MaskCalcVec::zero();
				case 1: return _mm256_set_epi32(0, 0, 0, 0, 0, 0, ~0, ~0);
				case 2: return _mm256_set_epi32(0, 0, 0, 0, ~0, ~0, ~0, ~0);
				default: return _mm256_set_epi32(0, 0, ~0, ~0, ~0, ~0, ~0, ~0);
			}
		}
	#endif
#elif VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER

	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		static vcp_inline MaskCalcVec vcp_simd_getInitMask(const size_t& i){
			static const MaskCalcVec possibleInitJMasks[VCP_VEC_SIZE] = { 0xFFFF, 0xFFFE, 0xFFFC, 0xFFF8, 0xFFF0, 0xFFE0, 0xFFC0, 0xFF80,
																		   0xFF00, 0xFE00, 0xFC00, 0xF800, 0xF000, 0xE000, 0xC000, 0x8000 };
			return possibleInitJMasks[i & static_cast<size_t>(VCP_VEC_SIZE_M1)];
		}

		static vcp_inline MaskCalcVec vcp_simd_getRemainderMask(const size_t& size) {
			static const MaskCalcVec possibleRemainderJMasks[VCP_VEC_SIZE] = { 0x0000, 0x0001, 0x0003, 0x0007, 0x000F, 0x001F, 0x003F, 0x007F,
																				0x00FF, 0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF };
			return possibleRemainderJMasks[size & static_cast<size_t>(VCP_VEC_SIZE_M1)];
		}
	#else /* VCP_DPDP */
		static vcp_inline MaskCalcVec vcp_simd_getInitMask(const size_t& i){
			static const MaskCalcVec possibleInitJMasks[VCP_VEC_SIZE] = { 0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE0, 0xC0, 0x80 };
			return possibleInitJMasks[i & static_cast<size_t>(VCP_VEC_SIZE_M1)];
		}

		static vcp_inline MaskCalcVec vcp_simd_getRemainderMask(const size_t& size) {
			static const MaskCalcVec possibleRemainderJMasks[VCP_VEC_SIZE] = { 0x00, 0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3F, 0x7F };
			return possibleRemainderJMasks[size & static_cast<size_t>(VCP_VEC_SIZE_M1)];
		}
	#endif
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
