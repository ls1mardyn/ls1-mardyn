/*
 * DoubleVec.h
 *
 *  Created on: 10 Nov 2016
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_DOUBLEVEC_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_DOUBLEVEC_H_

#include "./SIMD_TYPES.h"
#include "MaskVec.h"

namespace vcp {

class DoubleVec {
private:
	vcp_double_vec _d;

public:
	DoubleVec() {}

	operator vcp_double_vec() const {
		return _d;
	}

	DoubleVec(const vcp_double_vec & d) {
		_d = d;
	}

	static DoubleVec cast_MaskVec_to_DoubleVec(const MaskVec& m) {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<vcp_double_vec>(m);
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castsi128_pd(m);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castsi256_pd(m);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_castsi512_pd(m);
#endif
	}

	static MaskVec cast_DoubleVec_to_MaskVec(const DoubleVec& d) {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<vcp_mask_vec>(d);
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castpd_si128(d);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castpd_si256(d);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_castpd_si512(d);
#endif
	}

	static DoubleVec zero() {
		#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return 0.0;
		#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			return _mm_setzero_pd();
		#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			return _mm256_setzero_pd();
		#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return _mm512_setzero_pd();
		#endif
	}

	static DoubleVec ones() {
		return cast_MaskVec_to_DoubleVec(MaskVec::ones());
	}

	DoubleVec operator + (const DoubleVec& rhs) const {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d + rhs;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_add_pd(_d, rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_add_pd(_d, rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_add_pd(_d, rhs);
#endif
	}

	DoubleVec operator - (const DoubleVec& rhs) const {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d - rhs;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sub_pd(_d, rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sub_pd(_d, rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sub_pd(_d, rhs);
#endif
	}

	DoubleVec operator * (const DoubleVec& rhs) const {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d * rhs;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_mul_pd(_d, rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_mul_pd(_d, rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mul_pd(_d, rhs);
#endif
	}

	static DoubleVec fmadd(const DoubleVec & a, const DoubleVec& b, const DoubleVec& c ) {
#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) + c;
#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmadd_pd(a, b, c);
#else /* KNC or KNL */
		return _mm512_fmadd_pd(a, b, c);
#endif
	}

	static DoubleVec fnmadd(const DoubleVec & a, const DoubleVec& b, const DoubleVec& c ) {
#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return c - (a * b);
#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fnmadd_pd(a, b, c);
#else /* KNC or KNL */
		return _mm512_fnmadd_pd(a, b, c);
#endif
	}

	static DoubleVec fmsub(const DoubleVec & a, const DoubleVec& b, const DoubleVec& c) {
#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) - c;
#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmsub_pd(a, b, c);
#else /* KNC or KNL */
		return _mm512_fmsub_pd(a, b, c);
#endif
	}

	DoubleVec operator / (const DoubleVec& rhs) const {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d / rhs;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_div_pd(_d, rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_div_pd(_d, rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_div_pd(_d, rhs);
#endif
	}

	static DoubleVec sqrt (const DoubleVec& rhs) {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return sqrt(rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sqrt_pd(rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sqrt_pd(rhs);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sqrt_pd(rhs);
#endif
	}

	static DoubleVec scal_prod(
			const DoubleVec& a1, const DoubleVec& a2, const DoubleVec& a3,
			const DoubleVec& b1, const DoubleVec& b2, const DoubleVec& b3) {
		return fmadd(a1, b1, fmadd(a2, b2, a3 * b3));
	}

	static DoubleVec set1(const double& v) {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return v;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_set1_pd(v);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_set1_pd(v);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_set1_pd(a);
#endif
	}

	static DoubleVec aligned_load(const double * const a) {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return *a;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_load_pd(a);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_load_pd(a);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_load_pd(a);
#endif
	}

	static DoubleVec broadcast(const double * const a) {
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
	}

	void aligned_store(double * location) const {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		*location = _d;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		_mm_store_pd(location, _d);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		_mm256_store_pd(location, _d);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		_mm512_store_pd(location, _d);
#endif
	}

	//TODO: static functions returning object, not void members

	static DoubleVec unpack_lo(const DoubleVec& a, const DoubleVec& b) {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpacklo_pd(a,b);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpacklo_pd(a,b);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
#endif
	}

	static DoubleVec unpack_hi(const DoubleVec& a, const DoubleVec& b) {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpackhi_pd(a,b);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpackhi_pd(a, b);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
#endif
	}

	MaskVec operator < (const DoubleVec & rhs) const {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d < rhs;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_DoubleVec_to_MaskVec(_mm_cmplt_pd(_d, rhs));
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_DoubleVec_to_MaskVec(_mm256_cmp_pd(_d, rhs, _CMP_LT_OS));
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return cast_DoubleVec_to_MaskVec(_mm512_cmp_pd(_d, rhs, _CMP_LT_OS));
#endif
	}

	MaskVec operator != (const DoubleVec & rhs) const {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d != rhs;
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_DoubleVec_to_MaskVec(_mm_cmpeq_pd(_d, rhs));
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_DoubleVec_to_MaskVec(_mm256_cmp_pd(_d, rhs, _CMP_NEQ_OS));
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return cast_DoubleVec_to_MaskVec(_mm512_cmp_pd(_d, rhs, _CMP_NEQ_OS));
#endif
	}

	static DoubleVec apply_mask(const DoubleVec& d, const MaskVec& m) {
		return cast_MaskVec_to_DoubleVec(cast_DoubleVec_to_MaskVec(d) and m);
	}

	static DoubleVec aligned_load_mask(const double * const a, MaskVec m) {
#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return apply_mask(DoubleVec::aligned_load(a),m);
#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return apply_mask(DoubleVec::aligned_load(a),m);
#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_maskload_pd(a, m);
#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_load_pd(DoubleVec::zero(), m, a);
#endif
	}
};

} /* namespace vcp */

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_DOUBLEVEC_H_ */
