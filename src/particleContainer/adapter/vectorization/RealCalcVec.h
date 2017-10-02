/*
 * RealCalcVec.h
 *
 *  Created on: 10 Nov 2016
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALCALCVEC_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALCALCVEC_H_

#include "SIMD_TYPES.h"
#include "MaskVec.h"
#include "utils/mardyn_assert.h"
#include <fstream>
#include <cmath>

#include "utils/Logger.h"

namespace vcp {

class RealCalcVec {
private:
	vcp_real_calc_vec _d;

public:
	RealCalcVec() {}

	operator vcp_real_calc_vec() const {
		return _d;
	}

	RealCalcVec(const vcp_real_calc_vec & d) {
		_d = d;
	}

	static RealCalcVec cast_MaskVec_to_RealCalcVec(const MaskVec& m) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<vcp_real_calc_vec>(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castsi128_ps(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castsi256_ps(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0f / 0.0f); // do not use
	#endif
#else // VCP_DPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<vcp_real_calc_vec>(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castsi128_pd(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castsi256_pd(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0 / 0.0); // do not use
	#endif
#endif
	}

	static MaskVec cast_RealCalcVec_to_MaskVec(const RealCalcVec& d) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<vcp_mask_vec>(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castps_si128(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castps_si256(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return false; // do not use
	#endif
#else // VCP_DPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<vcp_mask_vec>(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castpd_si128(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castpd_si256(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return false; // do not use
	#endif
#endif
	}

	static RealCalcVec zero() {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0f;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_setzero_ps();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_setzero_ps();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_setzero_ps();
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_setzero_pd();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_setzero_pd();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_setzero_pd();
	#endif
#endif
	}

	static RealCalcVec ones() {
#if VCP_VEC_WIDTH != VCP_VEC_W_512
		return cast_MaskVec_to_RealCalcVec(MaskVec::ones());
#else
	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		return _mm512_castsi512_ps( _mm512_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0) );
	#else /* VCP_DPDP*/
		return _mm512_castsi512_pd( _mm512_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0) );
	#endif
#endif
	}

	RealCalcVec operator + (const RealCalcVec& rhs) const {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d + rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_add_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_add_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_add_ps(_d, rhs);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d + rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_add_pd(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_add_pd(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_add_pd(_d, rhs);
	#endif
#endif
	}

	RealCalcVec operator - (const RealCalcVec& rhs) const {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d - rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sub_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sub_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sub_ps(_d, rhs);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d - rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sub_pd(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sub_pd(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sub_pd(_d, rhs);
	#endif
#endif
	}

	RealCalcVec operator * (const RealCalcVec& rhs) const {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d * rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_mul_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_mul_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mul_ps(_d, rhs);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d * rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_mul_pd(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_mul_pd(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mul_pd(_d, rhs);
	#endif
#endif
	}

	static RealCalcVec fmadd(const RealCalcVec & a, const RealCalcVec& b, const RealCalcVec& c ) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) + c;
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmadd_ps(a, b, c);
	#else /* KNC or KNL */
		return _mm512_fmadd_ps(a, b, c);
	#endif
#else /* VCP_DPDP */
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) + c;
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmadd_pd(a, b, c);
	#else /* KNC or KNL */
		return _mm512_fmadd_pd(a, b, c);
	#endif
#endif
	}

	static RealCalcVec fnmadd(const RealCalcVec & a, const RealCalcVec& b, const RealCalcVec& c ) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return c - (a * b);
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fnmadd_ps(a, b, c);
	#else /* KNC or KNL */
		return _mm512_fnmadd_ps(a, b, c);
	#endif
#else /* VCP_DPDP */
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return c - (a * b);
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fnmadd_pd(a, b, c);
	#else /* KNC or KNL */
		return _mm512_fnmadd_pd(a, b, c);
	#endif
#endif
	}

	static RealCalcVec fmsub(const RealCalcVec & a, const RealCalcVec& b, const RealCalcVec& c) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) - c;
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmsub_ps(a, b, c);
	#else /* KNC or KNL */
		return _mm512_fmsub_ps(a, b, c);
	#endif
#else /* VCP_DPDP */
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) - c;
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmsub_pd(a, b, c);
	#else /* KNC or KNL */
		return _mm512_fmsub_pd(a, b, c);
	#endif
#endif
	}

	RealCalcVec operator / (const RealCalcVec& rhs) const {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d / rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_div_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_div_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_div_ps(_d, rhs);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d / rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_div_pd(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_div_pd(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_div_pd(_d, rhs);
	#endif
#endif
	}

	static RealCalcVec sqrt (const RealCalcVec& rhs) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return std::sqrt(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sqrt_ps(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sqrt_ps(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sqrt_ps(rhs);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return std::sqrt(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sqrt_pd(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sqrt_pd(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sqrt_pd(rhs);
	#endif
#endif
	}

	static RealCalcVec scal_prod(
		const RealCalcVec& a1, const RealCalcVec& a2, const RealCalcVec& a3,
		const RealCalcVec& b1, const RealCalcVec& b2, const RealCalcVec& b3) {
		return fmadd(a1, b1, fmadd(a2, b2, a3 * b3));
	}

	static RealCalcVec set1(const vcp_real_calc& v) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return v;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_set1_ps(v);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_set1_ps(v);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_set1_ps(v);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return v;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_set1_pd(v);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_set1_pd(v);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_set1_pd(v);
	#endif
#endif
	}

	static RealCalcVec aligned_load(const vcp_real_calc * const a) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return *a;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_load_ps(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_load_ps(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_load_ps(a);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return *a;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_load_pd(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_load_pd(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_load_pd(a);
	#endif
#endif
	}

	static RealCalcVec broadcast(const vcp_real_calc * const a) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return *a;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_load_ps1(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_broadcast_ss(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		#if VCP_VECTYPE==VCP_VEC_KNC or VCP_VECTYPE==VCP_VEC_KNC_GATHER
			return _mm512_extload_ps(a, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
		#else
			return _mm512_set1_ps(*a);
		#endif
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return *a;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_loaddup_pd(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_broadcast_sd(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		#if VCP_VECTYPE==VCP_VEC_KNC or VCP_VECTYPE==VCP_VEC_KNC_GATHER
			return _mm512_extload_pd(a, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
		#else
			return _mm512_set1_pd(*a);
		#endif
	#endif
#endif
	}

	void aligned_store(vcp_real_calc * location) const {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		*location = _d;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		_mm_store_ps(location, _d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		_mm256_store_ps(location, _d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		_mm512_store_ps(location, _d);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		*location = _d;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		_mm_store_pd(location, _d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		_mm256_store_pd(location, _d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		_mm512_store_pd(location, _d);
	#endif
#endif
	}

	static RealCalcVec unpack_lo(const RealCalcVec& a, const RealCalcVec& b) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0f / 0.0f; // makes no sense
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpacklo_ps(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpacklo_ps(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0f / 0.0f); // not used!
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0 / 0.0; // makes no sense
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpacklo_pd(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpacklo_pd(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0 / 0.0); // not used!
	#endif
#endif
	}

	static RealCalcVec unpack_hi(const RealCalcVec& a, const RealCalcVec& b) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0f / 0.0f; // makes no sense
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpackhi_ps(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpackhi_ps(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0f / 0.0f); // not used!
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0 / 0.0; // makes no sense
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpackhi_pd(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpackhi_pd(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0 / 0.0); // not used!
	#endif
#endif
	}

	MaskVec operator < (const RealCalcVec & rhs) const {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d < rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_RealCalcVec_to_MaskVec(_mm_cmplt_ps(_d, rhs));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_RealCalcVec_to_MaskVec(_mm256_cmp_ps(_d, rhs, _CMP_LT_OS));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cmp_ps_mask(_d, rhs, _CMP_LT_OS);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d < rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_RealCalcVec_to_MaskVec(_mm_cmplt_pd(_d, rhs));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_RealCalcVec_to_MaskVec(_mm256_cmp_pd(_d, rhs, _CMP_LT_OS));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cmp_pd_mask(_d, rhs, _CMP_LT_OS);
	#endif
#endif
	}

	MaskVec operator != (const RealCalcVec & rhs) const {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d != rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_RealCalcVec_to_MaskVec(_mm_cmpneq_ps(_d, rhs));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_RealCalcVec_to_MaskVec(_mm256_cmp_ps(_d, rhs, _CMP_NEQ_OS));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cmp_ps_mask(_d, rhs, _CMP_NEQ_UQ);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d != rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_RealCalcVec_to_MaskVec(_mm_cmpneq_pd(_d, rhs));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_RealCalcVec_to_MaskVec(_mm256_cmp_pd(_d, rhs, _CMP_NEQ_OS));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cmp_pd_mask(_d, rhs, _CMP_NEQ_UQ);
	#endif
#endif
	}

	/**
	 * \brief	Calculates 1/d and applies mask m on the result
	 * \detail	This could also be done by directly using operator / and method apply_mask(..),
	 * 			but this method is faster and more precise in some cases, because it calculates additional iterations of the Newton-Raphson method if beneficial.
	 * 			Furthermore, in DEBUG-mode it does an additional check for NaNs before returning the result.
	 */
	// Newton-Raphson division:
	// take the _mm*_rcp version with highest bit-precision
	// and iterate until the number of bits is >53 (DPDP) or >22 (SPSP, SPDP)
	//
	// An iteration of Newton-Raphson looks like this:
	// approximation * fnmadd(original_value, approximation, set1(2.0))
	//
	// AVX2 - some speed-up for single precision. For double prec slow down. Presumably because there are no dedicated fast intrinsics for packed double
	// AVX - slow-down: fma is not available; even more pressure on multiply-add imbalance? Leaving it out
	// KNC - small speed-up, because _mm512_div actually already uses Newton-Raphson, but doesn't assume that conversions double -> float are safe?
	// KNL - an educated guess would assume that AVX512ER is there for a reason :)
	static RealCalcVec fastReciprocal_mask(const RealCalcVec& d, const MaskVec& m) {

		/* reciprocal */
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const vcp_real_calc_vec inv_unmasked = _mm256_rcp_ps(d); //12bit
			#elif VCP_VEC_WIDTH == VCP_VEC_W_512
				#if VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
					const RealCalcVec inv_prec = _mm512_mask_rcp23_ps(zero(), m, d); //23bit
				#else /* VCP_VEC_KNL or VCP_VEC_KNL_GATHER */
					const RealCalcVec inv_prec = _mm512_maskz_rcp28_ps(m, d); //28bit
				#endif
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealCalcVec inv_unmasked = set1(1.0) / d;
			#endif
		#else /* VCP_DPDP */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const __m128 denom_ps = _mm256_cvtpd_ps(d);
				const __m128 recip_ps = _mm_rcp_ps(denom_ps);

				const vcp_real_calc_vec inv_unmasked = _mm256_cvtps_pd(recip_ps); //12bit
			#elif VCP_VEC_WIDTH == VCP_VEC_W_512
				#if VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
					// convert to single and truncate; no exceptions
				    const __m512 denom_ps = _mm512_mask_cvt_roundpd_pslo(_mm512_setzero_ps(), m, d, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
				    const __m512 D_recip_ps = _mm512_mask_rcp23_ps(_mm512_setzero_ps(), static_cast<__mmask16>(m), denom_ps);

				    const RealCalcVec inv = _mm512_mask_cvtpslo_pd(zero(), m, D_recip_ps); //23 bit
				#else /* VCP_VEC_KNL or VCP_VEC_KNL_GATHER */
					const RealCalcVec inv = _mm512_maskz_rcp28_pd(m, d); //28bit
				#endif
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealCalcVec inv_unmasked = set1(1.0) / d;
			#endif
		#endif

		/* mask and/or N-R-iterations if necessary */
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const RealCalcVec inv = apply_mask(inv_unmasked, m); //12bit

				const RealCalcVec inv_prec = inv * fnmadd(d, inv, set1(2.0)); //24bit, 1 N-R-Iteration
			#elif VCP_VEC_WIDTH == VCP_VEC_W_512
				//do nothing
				//no Newton-Raphson required in SP, as the rcp-intrinsic is already precise enough
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealCalcVec inv_prec = apply_mask(inv_unmasked, m);
			#endif
		#else /* VCP_DPDP */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const RealCalcVec inv = apply_mask(inv_unmasked, m); //12bit

				const RealCalcVec inv_24bits = inv * fnmadd(d, inv, set1(2.0)); 				//24bit, 1. N-R-Iteration
				const RealCalcVec inv_48bits = inv_24bits * fnmadd(d, inv_24bits, set1(2.0));	//48bit, 2. N-R-Iteration
				const RealCalcVec inv_prec = inv_48bits * fnmadd(d, inv_48bits, set1(2.0)); 	//96bit, 3. N-R-Iteration
			#elif VCP_VEC_WIDTH == VCP_VEC_W_512
				#if VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
				    const RealCalcVec inv_46bits = inv * fnmadd(d, inv, set1(2.0));				//46bit, 1. N-R-Iteration
				    const RealCalcVec inv_prec = inv_46bits * fnmadd(d, inv_46bits, set1(2.0)); //92bit, 2. N-R-Iteration
				#else /* VCP_VEC_KNL or VCP_VEC_KNL_GATHER */
				    const RealCalcVec inv_prec = inv * fnmadd(d, inv, set1(2.0)); //56bit, 1 N-R-Iteration
				#endif
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealCalcVec inv_prec = apply_mask(inv_unmasked, m);
			#endif
		#endif

		/* check for NaNs */
		#ifndef NDEBUG
			//set nan_found to zero only if NO NaN is found
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
				#if   VCP_VEC_WIDTH == VCP_VEC_W__64
					int nan_found = std::isnan(inv_prec) ? 1 : 0;
				#elif VCP_VEC_WIDTH == VCP_VEC_W_128
					//movemask-intrinsic is used to get a comparable int value out of vcp_real_calc_vec
					int nan_found = _mm_movemask_ps(_mm_cmpunord_ps(inv_prec, inv_prec));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_256
					__m128 low = _mm256_castps256_ps128(inv_prec);
					__m128 high = _mm256_extractf128_ps(inv_prec, 1);
					int nan_found = _mm_movemask_ps(_mm_cmpunord_ps(low, high));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_512
					//GCC doesnt know the convenient instruction, so the explicit comparison is used (at the moment)
					//vcp_mask_vec nan_found = _mm512_cmpunord_ps_mask(inv_prec, inv_prec);
					vcp_mask_vec nan_found = _mm512_cmp_ps_mask(inv_prec, inv_prec, _CMP_UNORD_Q);
				#endif
			#else /* VCP_DPDP */
				#if   VCP_VEC_WIDTH == VCP_VEC_W__64
					int nan_found = std::isnan(inv_prec) ? 1 : 0;
				#elif VCP_VEC_WIDTH == VCP_VEC_W_128
					int nan_found = _mm_movemask_pd(_mm_cmpunord_pd(inv_prec, inv_prec));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_256
					__m128d low = _mm256_castpd256_pd128(inv_prec);
					__m128d high = _mm256_extractf128_pd(inv_prec, 1);
					int nan_found = _mm_movemask_pd(_mm_cmpunord_pd(low, high));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_512
					vcp_mask_vec nan_found = _mm512_cmp_pd_mask(inv_prec, inv_prec, _CMP_UNORD_Q);
				#endif
			#endif
			if (nan_found) {
				Log::global_log->error() << "NaN detected! Perhaps forceMask is wrong" << std::endl;
				mardyn_assert(false);
			}
		#endif

		return inv_prec;
	} //fastReciprocal_mask(..)

	/**
	 * \brief	Calculates 1/sqrt(d) and applies mask m on the result
	 * \detail	This could also be done by directly using operator / and methods apply_mask(..), sqrt()
	 * 			but this method is faster because it makes use of special intrinsics.
	 * 			Makes use of Newton-Raphson iterations if necessary to acquire the necessary precision
	 * 			Furthermore, in DEBUG-mode it does an additional check for NaNs before returning the result.
	 */
	// Newton-Raphson division:
	// take the _mm*_rsqrt version with highest bit-precision
	// and iterate until the number of bits is >53 (DPDP) or >22 (SPSP, SPDP)
	//
	// An iteration of Newton-Raphson looks like this:
	// y(n+1) = y(n) * (1.5 - (d / 2) * y(n)²)	or
	// approximation * fnmadd(original_value * 0.5, approximation * approximation, set1(1.5))
	//
	// AVX2 - AVX2 - some speed-up for single precision. For double prec slow down. Presumably because there are no dedicated fast intrinsics for packed double
	// AVX - not supported
	// KNC - small speed-up, because _mm512_div actually already uses Newton-Raphson, but doesn't assume that conversions double -> float are safe?
	// KNL - an educated guess would assume that AVX512ER is there for a reason :)
	static RealCalcVec fastReciprocSqrt_mask(const RealCalcVec& d, const MaskVec& m) {

		/* reciprocal sqrt */
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const vcp_real_calc_vec invSqrt_unmasked = _mm256_rsqrt_ps(d); //12bit
			#elif VCP_VEC_WIDTH == VCP_VEC_W_512
				#if VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
					const RealCalcVec invSqrt_prec = _mm512_mask_rsqrt23_ps(zero(), m, d); //23bit
				#else /* VCP_VEC_KNL or VCP_VEC_KNL_GATHER */
					const RealCalcVec invSqrt_prec = _mm512_maskz_rsqrt28_ps(m, d); //28bit
				#endif
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealCalcVec invSqrt_unmasked = set1(1.0) / sqrt(d);
			#endif
		#else /* VCP_DPDP */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const __m128 denom_ps = _mm256_cvtpd_ps(d);
				const __m128 recipSqrt_ps = _mm_rsqrt_ps(denom_ps);

				const vcp_real_calc_vec invSqrt_unmasked = _mm256_cvtps_pd(recipSqrt_ps); //12bit
			#elif VCP_VEC_WIDTH == VCP_VEC_W_512
				#if VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
					// convert to single and truncate; no exceptions
				    const __m512 denom_ps = _mm512_mask_cvt_roundpd_pslo(_mm512_setzero_ps(), m, d, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
				    const __m512 D_recip_ps = _mm512_mask_rsqrt23_ps(_mm512_setzero_ps(), static_cast<__mmask16>(m), denom_ps);

				    const RealCalcVec invSqrt = _mm512_mask_cvtpslo_pd(zero(), m, D_recip_ps); //23 bit
				#else /* VCP_VEC_KNL or VCP_VEC_KNL_GATHER */
				    const RealCalcVec invSqrt = _mm512_maskz_rsqrt28_pd(m, d); //28bit
				#endif
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealCalcVec invSqrt_unmasked = set1(1.0) / sqrt(d);
			#endif
		#endif

		/* mask and/or N-R-iterations if necessary */
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const RealCalcVec invSqrt = apply_mask(invSqrt_unmasked, m); //12bit

				const RealCalcVec invSqrt_prec = invSqrt * fnmadd(d * set1(0.5), invSqrt * invSqrt, set1(1.5)); //24bit, 1 N-R-Iteration
			#elif VCP_VEC_WIDTH == VCP_VEC_W_512
				//do nothing
				//no Newton-Raphson required in SP, as the rsqrt-intrinsic is already precise enough
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealCalcVec invSqrt_prec = apply_mask(invSqrt_unmasked, m);
			#endif
		#else /* VCP_DPDP */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const RealCalcVec invSqrt = apply_mask(invSqrt_unmasked, m); //12bit

				const RealCalcVec d2 = d * set1(0.5);
				const RealCalcVec invSqrt_24bits = invSqrt 		  * fnmadd(d2, invSqrt * invSqrt, set1(1.5)); 				//24bit, 1. N-R-Iteration
				const RealCalcVec invSqrt_48bits = invSqrt_24bits * fnmadd(d2, invSqrt_24bits * invSqrt_24bits, set1(1.5));	//48bit, 2. N-R-Iteration
				const RealCalcVec invSqrt_prec   = invSqrt_48bits * fnmadd(d2, invSqrt_48bits * invSqrt_48bits, set1(1.5)); //96bit, 3. N-R-Iteration
			#elif VCP_VEC_WIDTH == VCP_VEC_W_512
				#if VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
					const RealCalcVec d2 = d * set1(0.5);
				    const RealCalcVec invSqrt_46bits = invSqrt * fnmadd(d2, invSqrt * invSqrt, set1(1.5));					  //46bit, 1. N-R-Iteration
				    const RealCalcVec invSqrt_prec = invSqrt_46bits * fnmadd(d2, invSqrt_46bits * invSqrt_46bits, set1(1.5)); //92bit, 2. N-R-Iteration
				#else /* VCP_VEC_KNL or VCP_VEC_KNL_GATHER */
				    const RealCalcVec invSqrt_prec = invSqrt * fnmadd(d * set1(0.5), invSqrt * invSqrt, set1(1.5)); //56bit, 1 N-R-Iteration
				#endif
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealCalcVec invSqrt_prec = apply_mask(invSqrt_unmasked, m);
			#endif
		#endif

		/* check for NaNs */
		#ifndef NDEBUG
			//set nan_found to zero only if NO NaN is found
			#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
				#if   VCP_VEC_WIDTH == VCP_VEC_W__64
					int nan_found = std::isnan(invSqrt_prec) ? 1 : 0;
				#elif VCP_VEC_WIDTH == VCP_VEC_W_128
					//movemask-intrinsic is used to get a comparable int value out of vcp_real_calc_vec
					int nan_found = _mm_movemask_ps(_mm_cmpunord_ps(invSqrt_prec, invSqrt_prec));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_256
					__m128 low = _mm256_castps256_ps128(invSqrt_prec);
					__m128 high = _mm256_extractf128_ps(invSqrt_prec, 1);
					int nan_found = _mm_movemask_ps(_mm_cmpunord_ps(low, high));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_512
					//GCC doesnt know the convenient instruction, so the explicit comparison is used (at the moment)
					//vcp_mask_vec nan_found = _mm512_cmpunord_ps_mask(invSqrt_prec, invSqrt_prec);
					vcp_mask_vec nan_found = _mm512_cmp_ps_mask(invSqrt_prec, invSqrt_prec, _CMP_UNORD_Q);
				#endif
			#else /* VCP_DPDP */
				#if   VCP_VEC_WIDTH == VCP_VEC_W__64
					int nan_found = std::isnan(invSqrt_prec) ? 1 : 0;
				#elif VCP_VEC_WIDTH == VCP_VEC_W_128
					int nan_found = _mm_movemask_pd(_mm_cmpunord_pd(invSqrt_prec, invSqrt_prec));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_256
					__m128d low = _mm256_castpd256_pd128(invSqrt_prec);
					__m128d high = _mm256_extractf128_pd(invSqrt_prec, 1);
					int nan_found = _mm_movemask_pd(_mm_cmpunord_pd(low, high));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_512
					vcp_mask_vec nan_found = _mm512_cmp_pd_mask(invSqrt_prec, invSqrt_prec, _CMP_UNORD_Q);
				#endif
			#endif
			if (nan_found) {
				Log::global_log->error() << "NaN detected! Perhaps forceMask is wrong" << std::endl;
				mardyn_assert(false);
			}
		#endif

		return invSqrt_prec;
	} //fastReciprocSqrt_mask(..)

	static RealCalcVec apply_mask(const RealCalcVec& d, const MaskVec& m) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if VCP_VEC_TYPE == VCP_NOVEC
		return m ? d : RealCalcVec::zero();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_mov_ps(RealCalcVec::zero(), m, d);
	#else // SSE, AVX, AVX2
		return cast_MaskVec_to_RealCalcVec(cast_RealCalcVec_to_MaskVec(d) and m);
	#endif
#else /* VCP_DPDP */
	#if VCP_VEC_TYPE == VCP_NOVEC
		return m ? d : RealCalcVec::zero();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_mov_pd(RealCalcVec::zero(), m, d);
	#else // SSE, AVX, AVX2
		return cast_MaskVec_to_RealCalcVec(cast_RealCalcVec_to_MaskVec(d) and m);
	#endif
#endif
	}

	static RealCalcVec aligned_load_mask(const vcp_real_calc * const a, MaskVec m) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return apply_mask(RealCalcVec::aligned_load(a),m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return apply_mask(RealCalcVec::aligned_load(a),m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_maskload_ps(a, m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_load_ps(RealCalcVec::zero(), m, a);
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return apply_mask(RealCalcVec::aligned_load(a),m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return apply_mask(RealCalcVec::aligned_load(a),m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_maskload_pd(a, m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_load_pd(RealCalcVec::zero(), m, a);
	#endif
#endif
	}

	vcp_inline
	static void horizontal_add_and_store(const RealCalcVec& a, vcp_real_calc * const mem_addr) {
#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		// add one float
		(*mem_addr) += a;

	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		// add four floats
		const RealCalcVec t1 = _mm_hadd_ps(a,a);
		const RealCalcVec t2 = _mm_hadd_ps(t1,t1);
		const vcp_real_calc t3 = _mm_cvtss_f32(t2);
		(*mem_addr) += t3;

	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		// add eight floats
		static const MaskVec memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 1<<31);
		const RealCalcVec t1 = _mm256_hadd_ps(a,a);
		const RealCalcVec t2 = _mm256_hadd_ps(t1,t1);
		const RealCalcVec t3 = _mm256_permute2f128_ps(t2, t2, 0x1);
		const RealCalcVec t4 = _mm256_add_ps(t2, t3);
		_mm256_maskstore_ps(
			mem_addr,
			memoryMask_first,
				t4 + RealCalcVec::aligned_load_mask(mem_addr, memoryMask_first)
		);


	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		// add sixteen floats
	#endif
#else /* VCP_DPDP */
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		// add one double
	    (*mem_addr) += a;

	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		// add two doubles
		const RealCalcVec t1 = _mm_hadd_pd(a, a);
		const RealCalcVec t2 = _mm_add_sd( t1, _mm_load_sd(mem_addr));
		_mm_store_sd( mem_addr, t2);

	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		// add four doubles
		static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
		const RealCalcVec a_t1 = _mm256_permute2f128_pd(a, a, 0x1);
		const RealCalcVec a_t2 = _mm256_hadd_pd(a, a_t1);
		const RealCalcVec a_t3 = _mm256_hadd_pd(a_t2, a_t2);
		_mm256_maskstore_pd(
			mem_addr,
			memoryMask_first,
				a_t3 + RealCalcVec::aligned_load_mask(mem_addr, memoryMask_first)
		);

	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		// add eight doubles
	#endif
#endif
	}
}; /* class RealCalcVec */

} /* namespace vcp */

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALCALCVEC_H_ */
