/*
 * RealVecFloat.h
 *
 *  Created on: 6 Feb 2018
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVECFLOAT_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVECFLOAT_H_

#include "RealVec.h"

// keep this file and RealVecDouble as close as possible, so that they can be examined via diff!

namespace vcp {

template<>
class RealVec<float> {
private:

	// own typedefs necessary
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		typedef float real_vec;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		typedef __m128 real_vec;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		typedef __m256 real_vec;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		typedef __m512 real_vec;
	#endif

protected:
	real_vec _d;

public:
	vcp_inline
	RealVec() {}

	vcp_inline
	operator real_vec() const {
		return _d;
	}

	vcp_inline
	RealVec(const real_vec & d) {
		_d = d;
	}

	vcp_inline
	RealVec(const RealVec& rhs) {
		_d = rhs._d;
	}

	vcp_inline
	RealVec& operator=(const RealVec& rhs) {
		_d = rhs._d;
		return *this;
	}

	vcp_inline
	static RealVec convertCalcToAccum(const RealVec & rcv) {
		return rcv;
	}

	vcp_inline
	static RealVec convertAccumToCalc(const RealVec & rav) {
		return rav;
	}

	vcp_inline
	static RealVec cvt_MaskVec_to_RealCalcVec(const MaskVec<float>& m) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<real_vec>(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_cvtepi32_ps(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_cvtepi32_ps(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cvtepi32_ps(m);
	#endif
	}

	vcp_inline
	static RealVec cast_MaskVec_to_RealCalcVec(const MaskVec<float>& m) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<real_vec>(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castsi128_ps(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castsi256_ps(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0f / 0.0f); // do not use
	#endif
	}

	vcp_inline
	static MaskVec<float> cast_RealCalcVec_to_MaskVec(const RealVec& d) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<MaskVec<float>>(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castps_si128(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castps_si256(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return false; // do not use
	#endif
	}

	vcp_inline
	static RealVec zero() {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0f;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_setzero_ps();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_setzero_ps();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_setzero_ps();
	#endif
	}

	vcp_inline
	static RealVec ones() {
#if VCP_VEC_WIDTH != VCP_VEC_W_512
		return cast_MaskVec_to_RealCalcVec(MaskVec<float>::ones());
#else
		return _mm512_castsi512_ps( _mm512_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0) );
#endif
	}

	vcp_inline
	RealVec operator + (const RealVec& rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d + rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_add_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_add_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_add_ps(_d, rhs);
	#endif
	}

	vcp_inline
	RealVec operator - (const RealVec& rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d - rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sub_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sub_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sub_ps(_d, rhs);
	#endif
	}

	vcp_inline
	RealVec operator * (const RealVec& rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d * rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_mul_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_mul_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mul_ps(_d, rhs);
	#endif
	}

	vcp_inline
	static RealVec fmadd(const RealVec & a, const RealVec& b, const RealVec& c ) {
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) + c;
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmadd_ps(a, b, c);
	#else /* KNL */
		return _mm512_fmadd_ps(a, b, c);
	#endif
	}

	vcp_inline
	static RealVec fnmadd(const RealVec & a, const RealVec& b, const RealVec& c ) {
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return c - (a * b);
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fnmadd_ps(a, b, c);
	#else /* KNL */
		return _mm512_fnmadd_ps(a, b, c);
	#endif
	}

	vcp_inline
	static RealVec fmsub(const RealVec & a, const RealVec& b, const RealVec& c) {
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) - c;
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmsub_ps(a, b, c);
	#else /* KNL */
		return _mm512_fmsub_ps(a, b, c);
	#endif
	}

	vcp_inline
	RealVec operator / (const RealVec& rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d / rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_div_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_div_ps(_d, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_div_ps(_d, rhs);
	#endif
	}

	vcp_inline
	static RealVec max (const RealVec& a, const RealVec& b) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return std::max(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_max_ps(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_max_ps(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_max_ps(a, b);
	#endif
	}

	vcp_inline
	static RealVec sqrt (const RealVec& rhs) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return std::sqrt(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sqrt_ps(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sqrt_ps(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sqrt_ps(rhs);
	#endif
	}

	vcp_inline
	static RealVec scal_prod(
		const RealVec& a1, const RealVec& a2, const RealVec& a3,
		const RealVec& b1, const RealVec& b2, const RealVec& b3) {
		return fmadd(a1, b1, fmadd(a2, b2, a3 * b3));
	}

	vcp_inline
	static RealVec set1(const float& v) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return v;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_set1_ps(v);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_set1_ps(v);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_set1_ps(v);
	#endif
	}

	vcp_inline
	static RealVec aligned_load(const float * const a) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return *a;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_load_ps(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_load_ps(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_load_ps(a);
	#endif
	}

	vcp_inline
	static RealVec broadcast(const float * const a) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return *a;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_load_ps1(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_broadcast_ss(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_set1_ps(*a);
	#endif
	}

	vcp_inline
	void aligned_store(float * location) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		*location = _d;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		_mm_store_ps(location, _d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		_mm256_store_ps(location, _d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		_mm512_store_ps(location, _d);
	#endif
	}

	vcp_inline
	static RealVec unpack_lo(const RealVec& a, const RealVec& b) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0f / 0.0f; // makes no sense
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpacklo_ps(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpacklo_ps(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0f / 0.0f); // not used!
	#endif
	}

	vcp_inline
	static RealVec unpack_hi(const RealVec& a, const RealVec& b) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0f / 0.0f; // makes no sense
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpackhi_ps(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpackhi_ps(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0f / 0.0f); // not used!
	#endif
	}

	vcp_inline
	MaskVec<float> operator < (const RealVec & rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d < rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_RealCalcVec_to_MaskVec(_mm_cmplt_ps(_d, rhs));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_RealCalcVec_to_MaskVec(_mm256_cmp_ps(_d, rhs, _CMP_LT_OS));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cmp_ps_mask(_d, rhs, _CMP_LT_OS);
	#endif
	}

	vcp_inline
	MaskVec<float> operator != (const RealVec & rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d != rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_RealCalcVec_to_MaskVec(_mm_cmpneq_ps(_d, rhs));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_RealCalcVec_to_MaskVec(_mm256_cmp_ps(_d, rhs, _CMP_NEQ_OS));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cmp_ps_mask(_d, rhs, _CMP_NEQ_UQ);
	#endif
	}

	vcp_inline
	static RealVec apply_mask(const RealVec& d, const MaskVec<float>& m) {
	#if VCP_VEC_TYPE == VCP_NOVEC
		return m ? d : RealVec::zero();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_mov_ps(RealVec::zero(), m, d);
	#else // SSE, AVX, AVX2
		return cast_MaskVec_to_RealCalcVec(cast_RealCalcVec_to_MaskVec(d) and m);
	#endif
	}

	vcp_inline
	static RealVec aligned_load_mask(const float * const a, MaskVec<float> m) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return apply_mask(RealVec::aligned_load(a),m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return apply_mask(RealVec::aligned_load(a),m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_maskload_ps(a, m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_load_ps(RealVec::zero(), m, a);
	#endif
	}

	vcp_inline
	static void horizontal_add_and_store(const RealVec& a, float * const mem_addr) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		// add one float
		(*mem_addr) += a;

	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		// add four floats
		const RealVec t1 = _mm_hadd_ps(a,a);
		const RealVec t2 = _mm_hadd_ps(t1,t1);
		const float t3 = _mm_cvtss_f32(t2);
		(*mem_addr) += t3;

	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		// add eight floats
		static const MaskVec<float> memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 1<<31);
		const RealVec t1 = _mm256_hadd_ps(a,a);
		const RealVec t2 = _mm256_hadd_ps(t1,t1);
		const RealVec t3 = _mm256_permute2f128_ps(t2, t2, 0x1);
		const RealVec t4 = _mm256_add_ps(t2, t3);
		_mm256_maskstore_ps(
			mem_addr,
			memoryMask_first,
				t4 + RealVec::aligned_load_mask(mem_addr, memoryMask_first)
		);


	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		// add sixteen floats
		//a		 = |  1 |  2 |  3 |  4 |  5 |....| 13 | 14 | 15 | 16 |
		//low	 = |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |
		//high	 = |  9 | 10 | 11 | 12 | 14 | 14 | 15 | 16 |
		__m256 low = _mm512_castps512_ps256(a);

		__m128 h1 = _mm512_extractf32x4_ps(a, 2);
		__m128 h2 = _mm512_extractf32x4_ps(a, 3);
		__m256 high = _mm256_castps128_ps256(h1);
		high = _mm256_insertf128_ps(high, h2, 1);

		//low+high 	= | 1+9 | 2+10 | 3+11 |...| 8+16 |
		//			= |  1  |  2   |  3   |...|  8   |
		//t1 = | 1+2 | 3+4 | 1+2 | 3+4 | 5+6 | 7+8 | 5+6 | 7+8 |
		__m256 t1 = _mm256_hadd_ps(low + high, low + high);
		//t1 = | 1234 | 1234 | 1234 | 1234 | 5678 | 5678 | 5678 | 5678 |
		t1 = _mm256_hadd_ps(t1, t1);
		//t2 = | 5678 | 5678 | 5678 | 5678 |
		__m128 t2 = _mm256_extractf128_ps(t1,1);
		//t3 = | 12345678 | 1234 | 1234 | 1234 |
		__m128 t3 = _mm_add_ss(_mm256_castps256_ps128(t1), t2);

		// only add first float value (12345678) of t3
		*mem_addr += _mm_cvtss_f32(t3);
	#endif
	}

	vcp_inline
	void aligned_load_add_store(float * location) const {
		RealVec dest = aligned_load(location);
		dest = dest + *this;
		dest.aligned_store(location);
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
	// KNL - an educated guess would assume that AVX512ER is there for a reason :)
	vcp_inline
	static RealVec fastReciprocal_mask(const RealVec& d, const MaskVec<float>& m) {

		/* reciprocal */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const real_vec inv_unmasked = _mm256_rcp_ps(d); //12bit
			#elif VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL
				const RealVec inv_prec = _mm512_maskz_rcp28_ps(m, d); //28bit
			#elif VCP_VEC_TYPE == VCP_VEC_AVX512F or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				const RealVec inv = _mm512_maskz_rcp14_ps(m, d); //14bit
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealVec inv_unmasked = set1(1.0) / d;
			#endif

		/* mask and/or N-R-iterations if necessary */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const RealVec inv = apply_mask(inv_unmasked, m); //12bit

				const RealVec inv_prec = inv * fnmadd(d, inv, set1(2.0)); //24bit, 1 N-R-Iteration
			#elif VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL
				//do nothing
				//no Newton-Raphson required in SP, as the rcp-intrinsic is already precise enough
			#elif VCP_VEC_TYPE == VCP_VEC_AVX512F or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				const RealVec inv_prec = inv * fnmadd(d, inv, set1(2.0)); // 28bit, 1 N-R-Iteration
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealVec inv_prec = apply_mask(inv_unmasked, m);
			#endif

		/* check for NaNs */
		#ifndef NDEBUG
			//set nan_found to zero only if NO NaN is found
				#if   VCP_VEC_WIDTH == VCP_VEC_W__64
					int nan_found = std::isnan(inv_prec) ? 1 : 0;
				#elif VCP_VEC_WIDTH == VCP_VEC_W_128
					//movemask-intrinsic is used to get a comparable int value out of real_vec
					int nan_found = _mm_movemask_ps(_mm_cmpunord_ps(inv_prec, inv_prec));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_256
					__m128 low = _mm256_castps256_ps128(inv_prec);
					__m128 high = _mm256_extractf128_ps(inv_prec, 1);
					int nan_found = _mm_movemask_ps(_mm_cmpunord_ps(low, high));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_512
					//GCC doesnt know the convenient instruction, so the explicit comparison is used (at the moment)
					//MaskVec<float> nan_found = _mm512_cmpunord_ps_mask(inv_prec, inv_prec);
					MaskVec<float> nan_found = _mm512_cmp_ps_mask(inv_prec, inv_prec, _CMP_UNORD_Q);
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
	// KNL - an educated guess would assume that AVX512ER is there for a reason :)
	vcp_inline
	static RealVec fastReciprocSqrt_mask(const RealVec& d, const MaskVec<float>& m) {

		/* reciprocal sqrt */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const real_vec invSqrt_unmasked = _mm256_rsqrt_ps(d); //12bit
			#elif VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL
				const RealVec invSqrt_prec = _mm512_maskz_rsqrt28_ps(m, d); //28bit
			#elif VCP_VEC_TYPE == VCP_VEC_AVX512F or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				const RealVec invSqrt = _mm512_maskz_rsqrt14_ps(m, d); //14bit
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealVec invSqrt_unmasked = set1(1.0) / sqrt(d);
			#endif

		/* mask and/or N-R-iterations if necessary */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const RealVec invSqrt = apply_mask(invSqrt_unmasked, m); //12bit

				const RealVec invSqrt_prec = invSqrt * fnmadd(d * set1(0.5), invSqrt * invSqrt, set1(1.5)); //24bit, 1 N-R-Iteration
			#elif VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL
				//do nothing
				//no Newton-Raphson required in SP, as the rsqrt-intrinsic is already precise enough
			#elif VCP_VEC_TYPE == VCP_VEC_AVX512F or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				const RealVec invSqrt_prec = invSqrt * fnmadd(d * set1(0.5), invSqrt * invSqrt, set1(1.5)); //28bit, 1 N-R-Iteration
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealVec invSqrt_prec = apply_mask(invSqrt_unmasked, m);
			#endif

		/* check for NaNs */
		#ifndef NDEBUG
			//set nan_found to zero only if NO NaN is found
				#if   VCP_VEC_WIDTH == VCP_VEC_W__64
					int nan_found = std::isnan(invSqrt_prec) ? 1 : 0;
				#elif VCP_VEC_WIDTH == VCP_VEC_W_128
					//movemask-intrinsic is used to get a comparable int value out of real_vec
					int nan_found = _mm_movemask_ps(_mm_cmpunord_ps(invSqrt_prec, invSqrt_prec));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_256
					__m128 low = _mm256_castps256_ps128(invSqrt_prec);
					__m128 high = _mm256_extractf128_ps(invSqrt_prec, 1);
					int nan_found = _mm_movemask_ps(_mm_cmpunord_ps(low, high));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_512
					//GCC doesnt know the convenient instruction, so the explicit comparison is used (at the moment)
					//MaskVec<float> nan_found = _mm512_cmpunord_ps_mask(invSqrt_prec, invSqrt_prec);
					MaskVec<float> nan_found = _mm512_cmp_ps_mask(invSqrt_prec, invSqrt_prec, _CMP_UNORD_Q);
				#endif
			if (nan_found) {
				Log::global_log->error() << "NaN detected! Perhaps forceMask is wrong" << std::endl;
				mardyn_assert(false);
			}
		#endif

		return invSqrt_prec;
	} //fastReciprocSqrt_mask(..)

}; /* class RealCalcVec */

} /* namespace vcp */

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVECFLOAT_H_ */
