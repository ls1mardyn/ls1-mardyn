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

#include "math.h"

/*
 * Check whether the file SIMD_TYPES.hpp has been included.
 * This file (SIMD_DEFINITIONS.hpp) needs some macros to be set properly by that file (SIMD_TYPES.hpp).
 */
#ifndef  SIMD_TYPES_H
    #error "SIMD_DEFINITIONS included without SIMD_TYPES! Never include this file directly! Include it only via SIMD_TYPES!"
#endif /* defined SIMD_TYPES_H */

#if VCP_VEC_TYPE==VCP_NOVEC
	static vcp_inline vcp_double_vec vcp_simd_zerov() { return 0.; }
	static vcp_inline vcp_double_vec vcp_simd_ones() { return 1.; }
	static const vcp_double_vec VCP_SIMD_ZEROV = vcp_simd_zerov();
	static const vcp_mask_vec VCP_SIMD_ZEROVM = false;
	static const vcp_mask_vec VCP_SIMD_ONESVM = true;
	static vcp_inline vcp_mask_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return a < b;}
	static vcp_inline vcp_mask_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return a == b;}
	static vcp_inline vcp_mask_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b) {return a != b;}

	/**
	 * do not use this to apply a mask, use vcp_simd_applymask instead !!!
	 * @param a
	 * @param b
	 * @return
	 */
	static vcp_inline vcp_mask_vec vcp_simd_and(const vcp_mask_vec& a, const vcp_mask_vec& b) {return a && b;}
	static vcp_inline vcp_mask_vec vcp_simd_or(const vcp_mask_vec& a, const vcp_mask_vec& b) {return a || b;}
	static vcp_inline vcp_mask_vec vcp_simd_xor(const vcp_mask_vec& a, const vcp_mask_vec& b) {return (a || b) && (not (a && b));}
	static vcp_inline vcp_double_vec vcp_simd_applymask(const vcp_double_vec& a, const vcp_mask_vec& mask) {return mask?a:0.;}

	static vcp_inline vcp_double_vec vcp_simd_add(const vcp_double_vec& a, const vcp_double_vec& b) {return a + b;}
	static vcp_inline vcp_double_vec vcp_simd_sub(const vcp_double_vec& a, const vcp_double_vec& b) {return a - b;}
	static vcp_inline vcp_double_vec vcp_simd_mul(const vcp_double_vec& a, const vcp_double_vec& b) {return a * b;}
	static vcp_inline vcp_double_vec vcp_simd_div(const vcp_double_vec& a, const vcp_double_vec& b) {return a / b;}
	static vcp_inline vcp_double_vec vcp_simd_sqrt(const vcp_double_vec& a) {return sqrt(a);}

	static vcp_inline vcp_double_vec vcp_simd_set1(const double& a) {return a;}
	static vcp_inline vcp_mask_vec vcp_simd_set1(const vcp_mask_single& a) {return a;}

	static vcp_inline vcp_double_vec vcp_simd_load(const double* const a) {return *a;}
	static vcp_inline vcp_mask_vec vcp_simd_load(const vcp_mask_single* const a) {return *a;}
	static vcp_inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return *a;}
	static vcp_inline void vcp_simd_store(double* location, const vcp_double_vec& a) {*location = a;}
	static vcp_inline void vcp_simd_store(vcp_mask_single* location, const vcp_mask_vec& a) {*location = a;}
	//static vcp_inline vcp_double_vec vcp_simd_unpacklo(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_unpacklo_pd(a,b);}
	//static vcp_inline vcp_double_vec vcp_simd_unpackhi(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_unpackhi_pd(a,b);}

	static vcp_inline bool vcp_simd_movemask(const vcp_mask_vec& a) {return a;}
	static vcp_inline vcp_double_vec vcp_simd_maskload(const double* const a, const vcp_mask_vec& mask) {return vcp_simd_applymask(vcp_simd_load(a),mask);}
	static vcp_inline vcp_mask_vec vcp_simd_getInitMask(const size_t& /*i*/){
		return true;
	}
	static vcp_inline vcp_mask_vec vcp_simd_getRemainderMask(const size_t& /*size*/){
		return false;
	}

#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	static vcp_inline vcp_double_vec vcp_simd_zerov() { return _mm_setzero_pd(); }
	static vcp_inline vcp_double_vec vcp_simd_ones() { return _mm_castsi128_pd( _mm_set_epi32(~0, ~0, ~0, ~0) ); }
	static const vcp_double_vec VCP_SIMD_ZEROV = vcp_simd_zerov();
	static const vcp_mask_vec VCP_SIMD_ZEROVM = _mm_castpd_si128(vcp_simd_zerov());
	static const vcp_mask_vec VCP_SIMD_ONESVM = _mm_castpd_si128(vcp_simd_ones());
	static vcp_inline vcp_mask_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_castpd_si128(_mm_cmplt_pd(a, b));}
	static vcp_inline vcp_mask_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_castpd_si128(_mm_cmpeq_pd(a, b));}
	static vcp_inline vcp_mask_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_castpd_si128(_mm_cmpneq_pd(a, b));}

	/**
	 * do not use this to apply a mask, use vcp_simd_applymask instead !!!
	 * @param a
	 * @param b
	 * @return
	 */
	static vcp_inline vcp_mask_vec vcp_simd_and(const vcp_mask_vec& a, const vcp_mask_vec& b) {return _mm_and_si128(a, b);}
	static vcp_inline vcp_mask_vec vcp_simd_or(const vcp_mask_vec& a, const vcp_mask_vec& b) {return _mm_or_si128(a, b);}
	static vcp_inline vcp_mask_vec vcp_simd_xor(const vcp_mask_vec& a, const vcp_mask_vec& b) {return _mm_xor_si128(a, b);}
	static vcp_inline vcp_double_vec vcp_simd_applymask(const vcp_double_vec& a, const vcp_mask_vec& mask) {return _mm_castsi128_pd(vcp_simd_and(_mm_castpd_si128(a), mask));}

	static vcp_inline vcp_double_vec vcp_simd_add(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_add_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_sub(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_sub_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_mul(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_mul_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_div(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_div_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_sqrt(const vcp_double_vec& a) {return _mm_sqrt_pd(a);}

	static vcp_inline vcp_double_vec vcp_simd_set1(const double& a) {return _mm_set1_pd(a);}
	static vcp_inline vcp_mask_vec vcp_simd_set1(const vcp_mask_single& a) {return _mm_set1_epi64x(a);}

	static vcp_inline vcp_double_vec vcp_simd_load(const double* const a) {return _mm_load_pd(a);}
	static vcp_inline vcp_mask_vec vcp_simd_load(const vcp_mask_single* const a) {return _mm_load_si128((const __m128i*)a);}
	static vcp_inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return _mm_loaddup_pd(a);}
	static vcp_inline void vcp_simd_store(double* location, const vcp_double_vec& a) {_mm_store_pd(location, a);}
	static vcp_inline void vcp_simd_store(vcp_mask_single* location, const vcp_mask_vec& a) {_mm_store_si128((__m128i*)location, a);}
	static vcp_inline vcp_double_vec vcp_simd_unpacklo(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_unpacklo_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_unpackhi(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_unpackhi_pd(a,b);}


	static vcp_inline bool vcp_simd_movemask(const vcp_mask_vec& a) {return _mm_movemask_epi8(a);}
	static vcp_inline vcp_double_vec vcp_simd_maskload(const double* const a, const vcp_mask_vec& mask) {return vcp_simd_applymask(vcp_simd_load(a),mask);}
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
	static vcp_inline vcp_double_vec vcp_simd_zerov() { return _mm256_setzero_pd(); }
	static vcp_inline vcp_double_vec vcp_simd_ones() { return _mm256_castsi256_pd( _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0) ); }
	static const vcp_double_vec VCP_SIMD_ZEROV = vcp_simd_zerov();
	static const vcp_mask_vec VCP_SIMD_ZEROVM = _mm256_castpd_si256(vcp_simd_zerov());
	static const vcp_mask_vec VCP_SIMD_ONESVM = _mm256_castpd_si256(vcp_simd_ones());
	static vcp_inline vcp_mask_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_LT_OS));}
	static vcp_inline vcp_mask_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_EQ_OS));}
	static vcp_inline vcp_mask_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_NEQ_OS));}
	/**
	 * do not use this to apply a mask, use vcp_simd_applymask instead !!!
	 * @param a
	 * @param b
	 * @return
	 */
	static vcp_inline vcp_mask_vec vcp_simd_and(const vcp_mask_vec& a, const vcp_mask_vec& b) {return _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(a), _mm256_castsi256_pd(b)));}
	static vcp_inline vcp_mask_vec vcp_simd_or(const vcp_mask_vec& a, const vcp_mask_vec& b) {return _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(a), _mm256_castsi256_pd(b)));}
	static vcp_inline vcp_mask_vec vcp_simd_xor(const vcp_mask_vec& a, const vcp_mask_vec& b) {return _mm256_castpd_si256(_mm256_xor_pd(_mm256_castsi256_pd(a), _mm256_castsi256_pd(b)));}
	static vcp_inline vcp_double_vec vcp_simd_applymask(const vcp_double_vec& a, const vcp_mask_vec& mask) {return _mm256_castsi256_pd(vcp_simd_and(_mm256_castpd_si256(a), mask));}

	static vcp_inline vcp_double_vec vcp_simd_add(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_add_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_sub(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_sub_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_mul(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_mul_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_div(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_div_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_sqrt(const vcp_double_vec& a) {return _mm256_sqrt_pd(a);}

	static vcp_inline vcp_double_vec vcp_simd_set1(const double& a) {return _mm256_set1_pd(a);}
	static vcp_inline vcp_mask_vec vcp_simd_set1(const vcp_mask_single& a) {return _mm256_set1_epi64x(a);}

	static vcp_inline vcp_double_vec vcp_simd_load(const double* const a) {return _mm256_load_pd(a);}
	static vcp_inline vcp_mask_vec vcp_simd_load(const vcp_mask_single* const a) {return _mm256_load_si256((const __m256i*)a);}
	static vcp_inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return _mm256_broadcast_sd(a);}
	static vcp_inline void vcp_simd_store(double* location, const vcp_double_vec& a) {_mm256_store_pd(location, a);}
	static vcp_inline void vcp_simd_store(vcp_mask_single* location, const vcp_mask_vec& a) {_mm256_store_si256((__m256i*)location, a);}
	static vcp_inline vcp_double_vec vcp_simd_unpacklo(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_unpacklo_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_unpackhi(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_unpackhi_pd(a,b);}

	static vcp_inline bool vcp_simd_movemask(const vcp_mask_vec& a) {return _mm256_movemask_pd(_mm256_castsi256_pd(a));}
	static vcp_inline vcp_double_vec vcp_simd_maskload(const double * const a, vcp_mask_vec mask) {return _mm256_maskload_pd(a, mask);}
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
			case 0: return VCP_SIMD_ZEROVM;
			case 1: return _mm256_set_epi32(0, 0, 0, 0, 0, 0, ~0, ~0);
			case 2: return _mm256_set_epi32(0, 0, 0, 0, ~0, ~0, ~0, ~0);
			default: return _mm256_set_epi32(0, 0, ~0, ~0, ~0, ~0, ~0, ~0);
		}
	}
#elif VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or\
	  VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER

	static vcp_inline vcp_double_vec vcp_simd_zerov() { return _mm512_castsi512_pd( _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ); }//exists
	static vcp_inline vcp_double_vec vcp_simd_ones() { return _mm512_castsi512_pd( _mm512_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0) ); }//exists
	static const vcp_double_vec VCP_SIMD_ZEROV = vcp_simd_zerov();
	static const vcp_mask_vec VCP_SIMD_ZEROVM = 0x00;
	static const vcp_mask_vec VCP_SIMD_ONESVM = 0xFF;
	//static vcp_inline vcp_double_vec mask_to_vcp_double_vec(const vcp_mask_vec& mask){return _mm512_mask_mov_pd(VCP_SIMD_ZEROV, mask, VCP_SIMD_ONESV);}
	static vcp_inline vcp_mask_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_cmp_pd_mask(a, b, _CMP_LT_OS);}//exists
	static vcp_inline vcp_mask_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_cmp_pd_mask(a, b, _CMP_EQ_OS);}//exists
	static vcp_inline vcp_mask_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b){return _mm512_cmp_pd_mask(a, b, _CMP_NEQ_UQ) ;}//exists
	/**
	 * do not use this to apply a mask, use vcp_simd_applymask instead !!!
	 * @param a
	 * @param b
	 * @return
	 */
	static vcp_inline vcp_mask_vec vcp_simd_and(const vcp_mask_vec& a, const vcp_mask_vec& b) {
		return a & b;
	}//only for mask vecs -> unsigned char or so...
	static vcp_inline vcp_mask_vec vcp_simd_or(const vcp_mask_vec& a, const vcp_mask_vec& b) {return a | b;}
	static vcp_inline vcp_mask_vec vcp_simd_xor(const vcp_mask_vec& a, const vcp_mask_vec& b) {return a ^ b;}
	static vcp_inline vcp_double_vec vcp_simd_applymask(const vcp_double_vec& a, const vcp_mask_vec& mask) {return _mm512_mask_mov_pd(VCP_SIMD_ZEROV, mask, a);}

	static vcp_inline vcp_double_vec vcp_simd_add(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_add_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_sub(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_sub_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_mul(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_mul_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_div(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_div_pd(a,b);}
	static vcp_inline vcp_double_vec vcp_simd_sqrt(const vcp_double_vec& a) {return _mm512_sqrt_pd(a);}

	static vcp_inline vcp_double_vec vcp_simd_set1(const double& a) {return _mm512_set1_pd(a);}//exists

	static vcp_inline vcp_double_vec vcp_simd_load(const double* const a) {return _mm512_load_pd(a);}
	static vcp_inline vcp_mask_vec vcp_simd_load(const vcp_mask_single* const a) {return *a;}
	static vcp_inline vcp_double_vec vcp_simd_maskload(const double * const a, vcp_mask_vec mask) {return _mm512_mask_load_pd(VCP_SIMD_ZEROV, mask, a);}
	//static vcp_inline vcp_lookupOrMask_vec vcp_simd_load(const vcp_lookupOrMask_single* const a) {return a;}

	#if VCP_VECTYPE==VCP_VEC_KNC or VCP_VECTYPE==VCP_VEC_KNC_GATHER
		static vcp_inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return _mm512_extload_pd(a, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);}
	#else
		static vcp_inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return _mm512_set1_pd(*a);}
	#endif

	static vcp_inline void vcp_simd_store(double* const location, const vcp_double_vec& a) {_mm512_store_pd(location, a);}
	static vcp_inline void vcp_simd_store(vcp_mask_single* const location, const vcp_mask_vec& a) {*location = a;}

	//static vcp_inline vcp_double_vec vcp_simd_unpacklo(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_unpacklo_pd(a,b);}//not needed
	//static vcp_inline vcp_double_vec vcp_simd_unpackhi(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_unpackhi_pd(a,b);}//not needed

	static vcp_inline bool vcp_simd_movemask(const vcp_mask_vec& a) {return a != VCP_SIMD_ZEROVM;}
	#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
		static vcp_inline vcp_lookupOrMask_vec vcp_simd_load(const vcp_lookupOrMask_single* const a) {return _mm512_load_epi64(a);}
		static vcp_inline void vcp_simd_store(vcp_lookupOrMask_single* location, const vcp_lookupOrMask_vec& a) {_mm512_store_epi64(location, a);}
	#endif

	static vcp_inline vcp_mask_vec vcp_simd_getInitMask(const size_t& i){
		static const vcp_mask_vec possibleInitJMasks[VCP_VEC_SIZE] = { 0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE0, 0xC0, 0x80 };
		return possibleInitJMasks[i & static_cast<size_t>(VCP_VEC_SIZE_M1)];
	}

	static vcp_inline vcp_mask_vec vcp_simd_getRemainderMask(const size_t& size) {
		static const vcp_mask_vec possibleRemainderJMasks[VCP_VEC_SIZE] = { 0x00, 0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3F, 0x7F };
		return possibleRemainderJMasks[size & static_cast<size_t>(VCP_VEC_SIZE_M1)];
	}
#endif

#if VCP_VEC_TYPE != VCP_NOVEC
	#ifdef __ICC
		#if __ICC < 1600 or VCP_VEC_TYPE == VCP_VEC_KNC or VCP_VEC_TYPE == VCP_VEC_KNC_GATHER
		//static vcp_inline vcp_double_vec operator < (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_lt(a, b); }//next three operators not compatible with gcc 4.7 or below
		//static vcp_inline vcp_double_vec operator == (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_eq(a, b); }
		//static vcp_inline vcp_double_vec operator != (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_neq(a, b); }
		//static vcp_inline vcp_double_vec operator & (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_and(a, b); } //the next three operators are not working with gcc at all
		//static vcp_inline vcp_double_vec operator | (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_or(a, b); }
		//static vcp_inline vcp_double_vec operator ^ (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_xor(a, b); }

		//#pragma message "icc commands for vectorization compiled."
        static vcp_inline vcp_double_vec operator + (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_add(a, b); }
        static vcp_inline vcp_double_vec operator - (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_sub(a, b); }
        static vcp_inline vcp_double_vec operator * (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_mul(a, b); }
        static vcp_inline vcp_double_vec operator / (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_div(a, b); }
		#else
		//#pragma message "icc commands skipped, since ICC >= 1600 and not MIC"
		#endif
	#else
		//#pragma message "icc commands for vectorization not compiled, since no icc detected."
	#endif
#else
	//#pragma message "icc commands for vectorization skipped, since novec"
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





// ------------- FMA, fmsub, gathers:

#if VCP_VEC_TYPE==VCP_VEC_AVX2
	static vcp_inline vcp_double_vec vcp_simd_fma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm256_fmadd_pd(a, b, c);
	}
	static vcp_inline vcp_double_vec vcp_simd_fnma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm256_fnmadd_pd(a, b, c);//-(a*b) + c
	}
	static vcp_inline vcp_double_vec vcp_simd_fms(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm256_fmsub_pd(a, b, c);
	}
#elif VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or \
	  VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
	static vcp_inline vcp_double_vec vcp_simd_fma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm512_fmadd_pd(a, b, c);
	}
	static vcp_inline vcp_double_vec vcp_simd_fnma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm512_fnmadd_pd(a, b, c);//-(a*b) + c
	}
	static vcp_inline vcp_double_vec vcp_simd_fms(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm512_fmsub_pd(a, b, c);
	}
#else //no fma available
	static vcp_inline vcp_double_vec vcp_simd_fma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return vcp_simd_add(vcp_simd_mul(a, b), c);
	}
	static vcp_inline vcp_double_vec vcp_simd_fnma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return vcp_simd_sub(c, vcp_simd_mul(a, b));//-(a*b) + c
	}
	static vcp_inline vcp_double_vec vcp_simd_fms(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return vcp_simd_sub(vcp_simd_mul(a, b), c);
	}
#endif

/**
 * calculates scalar product of a and b.
 * a1 * b1 + a2 * b2 + a3 * b3
 * @return
 */
static vcp_inline vcp_double_vec vcp_simd_scalProd(const vcp_double_vec& a1, const vcp_double_vec& a2, const vcp_double_vec& a3, const vcp_double_vec& b1, const vcp_double_vec& b2, const vcp_double_vec& b3) {
	return vcp_simd_fma(a1, b1, vcp_simd_fma(a2, b2, a3 * b3));
}














#endif /* SIMD_DEFINITIONS_H */
