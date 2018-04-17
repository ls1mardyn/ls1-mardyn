/*
 * RealAccumVec.h
 *
 *  Created on: 7 Feb 2018
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALACCUMVECBACKEND_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALACCUMVECBACKEND_H_

#include "RealVec.h"

namespace vcp {

#if VCP_VEC_WIDTH != VCP_VEC_W__64
// the novec case is handled differently, as it requires only one RealVec<double> to store its results.

class RealAccumVecBackend {

private:
	RealVec<double> _first;
	RealVec<double> _second;

public:
	vcp_inline
	RealAccumVecBackend() {}

	vcp_inline
	static RealAccumVecBackend convertCalcToAccum(const RealCalcVec & rcv) {
		RealVec<double> first = convert_low(rcv);
		RealVec<double> second = convert_high(rcv);
		return RealAccumVecBackend(first, second);
	}

	vcp_inline
	static RealCalcVec convertAccumToCalc(const RealAccumVecBackend & rav) {
		RealCalcVec ret = back_convert(rav._first, rav._second);
		return ret;
	}

	vcp_inline
	RealAccumVecBackend(const RealAccumVecBackend& rhs) {
		_first = rhs._first;
		_second = rhs._second;
	}

	vcp_inline
	RealAccumVecBackend(const RealVec<double>& first, const RealVec<double>& second) {
		_first = first;
		_second = second;
	}

	vcp_inline
	RealAccumVecBackend& operator=(const RealAccumVecBackend& rhs) {
		_first = rhs._first;
		_second = rhs._second;
		return *this;
	}

	vcp_inline
	static RealAccumVecBackend zero() {
		RealAccumVecBackend result;
		result._first = RealVec<double>::zero();
		result._second = RealVec<double>::zero();
		return result;
	}

	vcp_inline
	RealAccumVecBackend operator+(const RealAccumVecBackend& rhs) const {
		RealAccumVecBackend result;
		result._first = _first + rhs._first;
		result._second = _second + rhs._second;
		return result;
	}

	vcp_inline
	RealAccumVecBackend operator*(const RealAccumVecBackend& rhs) const {
		RealAccumVecBackend result;
		result._first = _first * rhs._first;
		result._second = _second * rhs._second;
		return result;
	}

	vcp_inline
	RealAccumVecBackend operator-(const RealAccumVecBackend& rhs) const {
		RealAccumVecBackend result;
		result._first = _first - rhs._first;
		result._second = _second - rhs._second;
		return result;
	}

	vcp_inline
	static RealAccumVecBackend fmadd(const RealAccumVecBackend & a, const RealAccumVecBackend& b, const RealAccumVecBackend& c ) {
		RealAccumVecBackend result;
		result._first = RealVec<double>::fmadd(a._first, b._first, c._first);
		result._second = RealVec<double>::fmadd(a._second, b._second, c._second);
		return result;
	}

	vcp_inline
	static RealAccumVecBackend fnmadd(const RealAccumVecBackend & a, const RealAccumVecBackend& b, const RealAccumVecBackend& c ) {
		RealAccumVecBackend result;
		result._first = RealVec<double>::fnmadd(a._first, b._first, c._first);
		result._second = RealVec<double>::fnmadd(a._second, b._second, c._second);
		return result;
	}

	vcp_inline
	void aligned_store(double * location) const {
		const size_t offset = sizeof(RealVec<double>) / sizeof(double);
		_first.aligned_store(location);
		_second.aligned_store(location + offset);
	}

	vcp_inline
	static RealAccumVecBackend aligned_load(const double * const a) {
		const size_t offset = sizeof(RealVec<double>) / sizeof(double);
		RealVec<double> first = RealVec<double>::aligned_load(a);
		RealVec<double> second = RealVec<double>::aligned_load(a + offset);
		return RealAccumVecBackend(first, second);
	}

	vcp_inline
	static RealAccumVecBackend aligned_load_mask(const double * const a, MaskVec<float> m) {
		// we need to make two masks of type MaskVec<double> from one MaskVec<float>
		MaskVec<double> m_lo, m_hi;
		convert_mask_vec(m, m_lo, m_hi);

		const size_t offset = sizeof(RealVec<double>) / sizeof(double);

		RealVec<double> first = RealVec<double>::aligned_load_mask(a, m_lo);
		RealVec<double> second = RealVec<double>::aligned_load_mask(a + offset, m_hi);

		return RealAccumVecBackend(first, second);
	}

	vcp_inline
	static RealAccumVecBackend set1(const double& v) {
		RealVec<double> first = RealVec<double>::set1(v);
		RealVec<double> second = RealVec<double>::set1(v);
		return RealAccumVecBackend(first, second);
	}

	vcp_inline
	static void horizontal_add_and_store(const RealAccumVecBackend& a, double * const mem_addr) {
		RealVec<double> sum = a._first + a._second;
		RealVec<double>::horizontal_add_and_store(sum, mem_addr);
	}

	vcp_inline
	void aligned_load_add_store(double * location) const {
		const size_t offset = sizeof(RealVec<double>) / sizeof(double);
		_first.aligned_load_add_store(location);
		_second.aligned_load_add_store(location + offset);
	}

	vcp_inline
	static RealAccumVecBackend scal_prod(
		const RealAccumVecBackend& a1, const RealAccumVecBackend& a2, const RealAccumVecBackend& a3,
		const RealAccumVecBackend& b1, const RealAccumVecBackend& b2, const RealAccumVecBackend& b3) {
		return fmadd(a1, b1, fmadd(a2, b2, a3 * b3));
	}


#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER
	vcp_inline
	static RealAccumVecBackend gather_load(const double * const src, const size_t& offset, const vcp_lookupOrMask_vec& lookup) {
		__m256i lookup_256i_lo = _mm512_extracti64x4_epi64(lookup, 0);
		__m256i lookup_256i_hi = _mm512_extracti64x4_epi64(lookup, 1);
		RealVec<double> first  (_mm512_i32gather_pd(lookup_256i_lo, src, 8));
		RealVec<double> second (_mm512_i32gather_pd(lookup_256i_hi, src, 8));
		return RealAccumVecBackend(first, second);
	}

	vcp_inline
	void gather_store(double* const addr, const size_t& offset, const vcp_lookupOrMask_vec& lookup) {
		__m256i lookup_256i_lo = _mm512_extracti64x4_epi64(lookup, 0);
		__m256i lookup_256i_hi = _mm512_extracti64x4_epi64(lookup, 1);
		_mm512_i32scatter_pd(addr, lookup_256i_lo, _first, 8);
		_mm512_i32scatter_pd(addr, lookup_256i_hi, _second, 8);
	}
#endif

	vcp_inline
	static RealVec<double> convert_low(const RealCalcVec& rhs) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		line not compiled
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_cvtps_pd(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_cvtps_pd(_mm256_extractf128_ps(rhs, 0));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(rhs), 0)));
	#endif
	}

	vcp_inline
	static RealVec<double> convert_high(const RealCalcVec& rhs) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		line not compiled
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(rhs), 8)));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_cvtps_pd(_mm256_extractf128_ps(rhs, 1));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(rhs), 1)));
	#endif
	}

	vcp_inline
	static RealCalcVec back_convert(const RealVec<double>& lo, const RealVec<double>& hi) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		line not compiled
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		__m128 c_lo = _mm_cvtpd_ps(lo);
		__m128 c_hi = _mm_cvtpd_ps(hi);
		return _mm_movelh_ps(c_lo, c_hi);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		__m128 c_lo = _mm256_cvtpd_ps(lo);
		__m128 c_hi = _mm256_cvtpd_ps(hi);

		__m256 ret = _mm256_castps128_ps256(c_lo);
		ret = _mm256_insertf128_ps(ret, c_hi, 1);

		return ret;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		__m256 c_lo = _mm512_cvtpd_ps(lo);
		__m256 c_hi = _mm512_cvtpd_ps(hi);

		__m512 ret = _mm512_castps256_ps512(c_lo);
		ret = _mm512_insertf32x8(ret, c_hi, 1);

		return ret;
	#endif
	}

	vcp_inline
	static void convert_mask_vec(const MaskVec<float>& src, MaskVec<double>& lo, MaskVec<double>& hi) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		line not compiled
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128

		lo = _mm_unpacklo_epi32(src, src);
		hi = _mm_unpackhi_epi32(src, src);

	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		__m256 v_3210 = _mm256_castsi256_ps(src);

		__m256i v_2200 = _mm256_castps_si256(_mm256_unpacklo_ps(v_3210, v_3210));
		__m256i v_3311 = _mm256_castps_si256(_mm256_unpackhi_ps(v_3210, v_3210));

		// need to swap 11 and 22
		auto v_10 = _mm256_extract_epi64 (v_3311, 0);
		auto v_11 = _mm256_extract_epi64 (v_3311, 1);

		auto v_20 = _mm256_extract_epi64 (v_2200, 2);
		auto v_21 = _mm256_extract_epi64 (v_2200, 3);

		__m256i v_2100 = _mm256_insert_epi64 (v_2200, v_10, 2);
		__m256i v_1100 = _mm256_insert_epi64 (v_2100, v_11, 3);
		lo = v_1100;

		__m256i v_3312 = _mm256_insert_epi64 (v_3311, v_20, 0);
		__m256i v_3322 = _mm256_insert_epi64 (v_3312, v_21, 1);
		hi = v_3322;

	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		// need to make two __mmask8 from one __mmask16
		// the intrinsics are not very helpful for working with mmask*, so...
	    union {
	        __mmask16 _wide;
	        __mmask8 _narrow[2];
	    } merged;
	    merged._wide = src;
	    lo = merged._narrow[0];
	    hi = merged._narrow[1];

	#endif
	}
};

#elif VCP_VEC_WIDTH == VCP_VEC_W__64

class RealAccumVecBackend : public RealVec<double> {
public:
	vcp_inline
	RealAccumVecBackend() {}

	vcp_inline
	RealAccumVecBackend(const RealVec<double>& rcv) : RealVec<double>() {
		this->_d = rcv;
	}

	vcp_inline
	static RealAccumVecBackend convertCalcToAccum(const RealCalcVec & rcv) {
		RealAccumVecBackend result;
		result._d = rcv;
		return result;
	}

	vcp_inline
	static RealCalcVec convertAccumToCalc(const RealAccumVecBackend & rav) {
		RealCalcVec result(rav);
		return result;
	}

	vcp_inline
	static RealAccumVecBackend aligned_load_mask(const double * const a, MaskVec<float> m) {
		return apply_mask(aligned_load(a),MaskVec<double>(m));
	}
};

#endif /* VCP_VEC_WIDTH */

} /* namespace vcp */

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALACCUMVECBACKEND_H_ */
