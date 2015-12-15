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

#include "math.h"

/*
 * Check whether the file SIMD_TYPES.hpp has been included.
 * This file (SIMD_DEFINITIONS.hpp) needs some macros to be set properly by that file (SIMD_TYPES.hpp).
 */
#ifndef  SIMD_TYPES_H
    #error "SIMD_DEFINITIONS included without SIMD_TYPES! Never include this file directly! Include it only via SIMD_TYPES!"
#endif /* defined SIMD_TYPES_H */

#if VCP_VEC_TYPE==VCP_NOVEC
	static inline vcp_double_vec vcp_simd_zerov() { return 0.; }
	static inline vcp_double_vec vcp_simd_ones() { return 1.; }
	static const vcp_double_vec VCP_SIMD_ZEROV = vcp_simd_zerov();
	static const vcp_double_vec VCP_SIMD_ONESV = vcp_simd_ones();
	static inline vcp_double_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return a < b;}
	static inline vcp_double_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return a == b;}
	static inline vcp_double_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b) {return a != b;}

	/**
	 * do not use this to apply a mask, use vcp_simd_applymask instead !!!
	 * @param a
	 * @param b
	 * @return
	 */
	static inline vcp_double_vec vcp_simd_and(const vcp_double_vec& a, const vcp_double_vec& b) {return a && b;}
	static inline vcp_double_vec vcp_simd_or(const vcp_double_vec& a, const vcp_double_vec& b) {return a || b;}
	static inline vcp_double_vec vcp_simd_xor(const vcp_double_vec& a, const vcp_double_vec& b) {return (a || b) && (not (a && b));}
	static inline vcp_double_vec vcp_simd_applymask(const vcp_double_vec& a, const vcp_double_vec& mask) {return mask?a:0l;}

	static inline vcp_double_vec vcp_simd_add(const vcp_double_vec& a, const vcp_double_vec& b) {return a + b;}
	static inline vcp_double_vec vcp_simd_sub(const vcp_double_vec& a, const vcp_double_vec& b) {return a - b;}
	static inline vcp_double_vec vcp_simd_mul(const vcp_double_vec& a, const vcp_double_vec& b) {return a * b;}
	static inline vcp_double_vec vcp_simd_div(const vcp_double_vec& a, const vcp_double_vec& b) {return a / b;}
	static inline vcp_double_vec vcp_simd_sqrt(const vcp_double_vec& a) {return sqrt(a);}

	static inline vcp_double_vec vcp_simd_set1(const double& a) {return a;}

	static inline vcp_double_vec vcp_simd_load(const double* const a) {return *a;}
	static inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return *a;}
	static inline void vcp_simd_store(double* location, const vcp_double_vec& a) {*location = a;}
	//static inline vcp_double_vec vcp_simd_unpacklo(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_unpacklo_pd(a,b);}
	//static inline vcp_double_vec vcp_simd_unpackhi(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_unpackhi_pd(a,b);}

	static inline int vcp_simd_movemask(const vcp_double_vec& a) {return a?1:0;}
	static inline int vcp_simd_movemask(const vcp_doublesizedmask_vec& a) {return a?1:0;}

#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	static inline vcp_double_vec vcp_simd_zerov() { return _mm_setzero_pd(); }
	static inline vcp_double_vec vcp_simd_ones() { return _mm_castsi128_pd( _mm_set_epi32(~0, ~0, ~0, ~0) ); }
	static const vcp_double_vec VCP_SIMD_ZEROV = vcp_simd_zerov();
	static const vcp_double_vec VCP_SIMD_ONESV = vcp_simd_ones();
	static inline vcp_double_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_cmplt_pd(a, b);}
	static inline vcp_double_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_cmpeq_pd(a, b);}
	static inline vcp_double_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_cmpneq_pd(a, b);}

	/**
	 * do not use this to apply a mask, use vcp_simd_applymask instead !!!
	 * @param a
	 * @param b
	 * @return
	 */
	static inline vcp_double_vec vcp_simd_and(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_and_pd(a, b);}
	static inline vcp_double_vec vcp_simd_or(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_or_pd(a, b);}
	static inline vcp_double_vec vcp_simd_xor(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_xor_pd(a, b);}
	static inline vcp_double_vec vcp_simd_applymask(const vcp_double_vec& a, const vcp_double_vec& mask) {return vcp_simd_and(a, mask);}

	static inline vcp_double_vec vcp_simd_add(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_add_pd(a,b);}
	static inline vcp_double_vec vcp_simd_sub(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_sub_pd(a,b);}
	static inline vcp_double_vec vcp_simd_mul(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_mul_pd(a,b);}
	static inline vcp_double_vec vcp_simd_div(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_div_pd(a,b);}
	static inline vcp_double_vec vcp_simd_sqrt(const vcp_double_vec& a) {return _mm_sqrt_pd(a);}

	static inline vcp_double_vec vcp_simd_set1(const double& a) {return _mm_set1_pd(a);}

	static inline vcp_double_vec vcp_simd_load(const double* const a) {return _mm_load_pd(a);}
	static inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return _mm_loaddup_pd(a);}
	static inline void vcp_simd_store(double* location, const vcp_double_vec& a) {_mm_store_pd(location, a);}
	static inline vcp_double_vec vcp_simd_unpacklo(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_unpacklo_pd(a,b);}
	static inline vcp_double_vec vcp_simd_unpackhi(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm_unpackhi_pd(a,b);}


	static inline int vcp_simd_movemask(const vcp_double_vec& a) {return _mm_movemask_pd(a);}

#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2
	static inline vcp_double_vec vcp_simd_zerov() { return _mm256_setzero_pd(); }
	static inline vcp_double_vec vcp_simd_ones() { return _mm256_castsi256_pd( _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0) ); }
	static const vcp_double_vec VCP_SIMD_ZEROV = vcp_simd_zerov();
	static const vcp_double_vec VCP_SIMD_ONESV = vcp_simd_ones();
	static inline vcp_double_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_cmp_pd(a, b, _CMP_LT_OS);}
	static inline vcp_double_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_cmp_pd(a, b, _CMP_EQ_OS);}
	static inline vcp_double_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_cmp_pd(a, b, _CMP_NEQ_OS);}
	/**
	 * do not use this to apply a mask, use vcp_simd_applymask instead !!!
	 * @param a
	 * @param b
	 * @return
	 */
	static inline vcp_double_vec vcp_simd_and(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_and_pd(a, b);}
	static inline vcp_double_vec vcp_simd_or(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_or_pd(a, b);}
	static inline vcp_double_vec vcp_simd_xor(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_xor_pd(a, b);}
	static inline vcp_double_vec vcp_simd_applymask(const vcp_double_vec& a, const vcp_double_vec& mask) {return vcp_simd_and(a, mask);}

	static inline vcp_double_vec vcp_simd_add(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_add_pd(a,b);}
	static inline vcp_double_vec vcp_simd_sub(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_sub_pd(a,b);}
	static inline vcp_double_vec vcp_simd_mul(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_mul_pd(a,b);}
	static inline vcp_double_vec vcp_simd_div(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_div_pd(a,b);}
	static inline vcp_double_vec vcp_simd_sqrt(const vcp_double_vec& a) {return _mm256_sqrt_pd(a);}

	static inline vcp_double_vec vcp_simd_set1(const double& a) {return _mm256_set1_pd(a);}

	static inline vcp_double_vec vcp_simd_load(const double* const a) {return _mm256_load_pd(a);}
	static inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return _mm256_broadcast_sd(a);}
	static inline void vcp_simd_store(double* location, const vcp_double_vec& a) {_mm256_store_pd(location, a);}
	static inline vcp_double_vec vcp_simd_unpacklo(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_unpacklo_pd(a,b);}
	static inline vcp_double_vec vcp_simd_unpackhi(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_unpackhi_pd(a,b);}

	static inline int vcp_simd_movemask(const vcp_double_vec& a) {return _mm256_movemask_pd(a);}
	static inline vcp_double_vec vcp_simd_maskload(double const * a, vcp_mask_vec b) {return _mm256_maskload_pd(a, b);}
#elif VCP_VEC_TYPE==VCP_VEC_MIC

	static inline vcp_double_vec vcp_simd_zerov() { return _mm512_castsi512_pd( _mm512_set_epi32(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) ); }//exists
	static inline vcp_double_vec vcp_simd_ones() { return _mm512_castsi512_pd( _mm512_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0) ); }//exists
	static const vcp_double_vec VCP_SIMD_ZEROV = vcp_simd_zerov();
	static const vcp_double_vec VCP_SIMD_ONESV = vcp_simd_ones();
	static inline vcp_double_vec mask_to_vcp_double_vec(const vcp_mask_vec& mask){return _mm512_mask_mov_pd(VCP_SIMD_ZEROV, mask, VCP_SIMD_ONESV);}
	static inline vcp_double_vec vcp_simd_lt(const vcp_double_vec& a, const vcp_double_vec& b) {return mask_to_vcp_double_vec(_mm512_cmp_pd_mask(a, b, _CMP_LT_OS));}//exists
	static inline vcp_double_vec vcp_simd_eq(const vcp_double_vec& a, const vcp_double_vec& b) {return mask_to_vcp_double_vec(_mm512_cmp_pd_mask(a, b, _CMP_EQ_OS));}//exists
	static inline vcp_double_vec vcp_simd_neq(const vcp_double_vec& a, const vcp_double_vec& b){
		return mask_to_vcp_double_vec(
			_mm512_cmp_pd_mask(
					a,
					b,
					_CMP_NEQ_UQ)
		);}//exists
	/**
	 * do not use this to apply a mask, use vcp_simd_applymask instead !!!
	 * @param a
	 * @param b
	 * @return
	 */
	static inline vcp_double_vec vcp_simd_and(const vcp_double_vec& a, const vcp_double_vec& b) {
		return _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(a), _mm512_castpd_si512(b)));
	}//only for mask vecs -> unsigned char or so...
	static inline vcp_double_vec vcp_simd_or(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_castsi512_pd(_mm512_or_epi64(_mm512_castpd_si512(a), _mm512_castpd_si512(b)));}
	static inline vcp_double_vec vcp_simd_xor(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_castsi512_pd(_mm512_xor_epi64(_mm512_castpd_si512(a), _mm512_castpd_si512(b)));}
	static inline vcp_double_vec vcp_simd_applymask(const vcp_double_vec& a, const vcp_double_vec& mask) {return vcp_simd_and(a, mask);}

	static inline vcp_double_vec vcp_simd_add(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_add_pd(a,b);}
	static inline vcp_double_vec vcp_simd_sub(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_sub_pd(a,b);}
	static inline vcp_double_vec vcp_simd_mul(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_mul_pd(a,b);}
	static inline vcp_double_vec vcp_simd_div(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm512_div_pd(a,b);}
	static inline vcp_double_vec vcp_simd_sqrt(const vcp_double_vec& a) {return _mm512_sqrt_pd(a);}

	static inline vcp_double_vec vcp_simd_set1(const double& a) {return _mm512_set1_pd(a);}//exists

	static inline vcp_double_vec vcp_simd_load(const double* const a) {return _mm512_load_pd(a);}
	static inline vcp_double_vec vcp_simd_broadcast(const double* const a) {return _mm512_extload_pd(a, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);}
	static inline void vcp_simd_store(double* location, const vcp_double_vec& a) {_mm512_store_pd(location, a);}
	//static inline vcp_double_vec vcp_simd_unpacklo(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_unpacklo_pd(a,b);}//not needed
	//static inline vcp_double_vec vcp_simd_unpackhi(const vcp_double_vec& a, const vcp_double_vec& b) {return _mm256_unpackhi_pd(a,b);}//not needed

	static inline int vcp_simd_movemask(const vcp_double_vec& a) {return (int) _mm512_cmp_pd_mask(a, VCP_SIMD_ZEROV, _CMP_NEQ_UQ);}
#endif

#if VCP_VEC_TYPE != VCP_NOVEC
	#ifdef __ICC
		//static inline vcp_double_vec operator < (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_lt(a, b); }//next three operators not compatible with gcc 4.7 or below
		//static inline vcp_double_vec operator == (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_eq(a, b); }
		//static inline vcp_double_vec operator != (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_neq(a, b); }
		//static inline vcp_double_vec operator & (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_and(a, b); } //the next three operators are not working with gcc at all
		//static inline vcp_double_vec operator | (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_or(a, b); }
		//static inline vcp_double_vec operator ^ (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_xor(a, b); }


        static inline vcp_double_vec operator + (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_add(a, b); }
        static inline vcp_double_vec operator - (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_sub(a, b); }
        static inline vcp_double_vec operator * (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_mul(a, b); }
        static inline vcp_double_vec operator / (const vcp_double_vec& a, const vcp_double_vec& b) { return vcp_simd_div(a, b); }
	#endif
#endif

/**
 * if num is not divisible by VCP_VEC_SIZE finds the next bigger multiple of VCP_VEC_SIZE, otherwise returns num.
 * @param num
 * @return
 */
template<class T>
static inline T vcp_ceil_to_vec_size(const T& num){
	return (num + static_cast<T>(VCP_VEC_SIZE_M1)) & (~static_cast<T>(VCP_VEC_SIZE_M1));
}

/**
 * if num is not divisible by VCP_VEC_SIZE finds the next smaller multiple of VCP_VEC_SIZE, otherwise returns num.
 * @param num
 * @return
 */
template<class T>
static inline T vcp_floor_to_vec_size(const T& num){
	return num & (~static_cast<T>(VCP_VEC_SIZE_M1));
}





// ------------- FMA, fmsub, gathers:

#if VCP_VEC_TYPE==VCP_VEC_AVX2
	static inline vcp_double_vec vcp_simd_fma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm256_fmadd_pd(a, b, c);
	}
	static inline vcp_double_vec vcp_simd_fnma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm256_fnmadd_pd(a, b, c);//-(a*b) + c
	}
	static inline vcp_double_vec vcp_simd_fms(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm256_fmsub_pd(a, b, c);
	}
#elif VCP_VEC_TYPE==VCP_VEC_MIC
	static inline vcp_double_vec vcp_simd_fma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm512_fmadd_pd(a, b, c);
	}
	static inline vcp_double_vec vcp_simd_fnma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm512_fnmadd_pd(a, b, c);//-(a*b) + c
	}
	static inline vcp_double_vec vcp_simd_fms(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return _mm512_fmsub_pd(a, b, c);
	}
#else //no fma available
	static inline vcp_double_vec vcp_simd_fma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return vcp_simd_add(vcp_simd_mul(a, b), c);
	}
	static inline vcp_double_vec vcp_simd_fnma(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return vcp_simd_sub(c, vcp_simd_mul(a, b));//-(a*b) + c
	}
	static inline vcp_double_vec vcp_simd_fms(const vcp_double_vec& a, const vcp_double_vec& b, const vcp_double_vec& c) {
		return vcp_simd_sub(vcp_simd_mul(a, b), c);
	}
#endif

/**
 * calculates scalar product of a and b.
 * a1 * b1 + a2 * b2 + a3 * b3
 * @return
 */
static inline vcp_double_vec vcp_simd_scalProd(const vcp_double_vec& a1, const vcp_double_vec& a2, const vcp_double_vec& a3, const vcp_double_vec& b1, const vcp_double_vec& b2, const vcp_double_vec& b3) {
	return vcp_simd_fma(a1, b1, vcp_simd_fma(a2, b2, a3 * b3));
}














#endif /* SIMD_DEFINITIONS_H */
