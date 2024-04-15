/*
 * RealVecDouble.h
 *
 *  Created on: 6 Feb 2018
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVECDOUBLE_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVECDOUBLE_H_

#include "RealVec.h"

// keep this file and RealVecFloat as close as possible, so that they can be examined via diff!

namespace vcp {

template<>
class RealVec<double> {
private:

	// own typedefs necessary
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		typedef double real_vec;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		typedef __m128d real_vec;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		typedef __m256d real_vec;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		typedef __m512d real_vec;
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
	static RealVec cvt_MaskVec_to_RealCalcVec(const MaskVec<double>& m) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<real_vec>(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_cvtepi64_pd(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_cvtepi64_pd(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cvtepi64_pd(m);
	#endif
	}

	vcp_inline
	static RealVec cast_MaskVec_to_RealCalcVec(const MaskVec<double>& m) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<real_vec>(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castsi128_pd(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castsi256_pd(m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0 / 0.0); // do not use
	#endif
	}

	vcp_inline
	static MaskVec<double> cast_RealCalcVec_to_MaskVec(const RealVec& d) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return static_cast<MaskVec<double>>(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_castpd_si128(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_castpd_si256(d);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return false; // do not use
	#endif
	}

	vcp_inline
	static RealVec zero() {
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

	vcp_inline
	static RealVec ones() {
#if VCP_VEC_WIDTH != VCP_VEC_W_512
		return cast_MaskVec_to_RealCalcVec(MaskVec<double>::ones());
#else
		return _mm512_castsi512_pd( _mm512_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0) );
#endif
	}

	vcp_inline
	RealVec operator + (const RealVec& rhs) const {
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

	vcp_inline
	RealVec operator - (const RealVec& rhs) const {
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

	vcp_inline
	RealVec operator * (const RealVec& rhs) const {
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

	vcp_inline
	static RealVec fmadd(const RealVec & a, const RealVec& b, const RealVec& c ) {
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) + c;
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmadd_pd(a, b, c);
	#else /* KNL */
		return _mm512_fmadd_pd(a, b, c);
	#endif
	}

	vcp_inline
	static RealVec fnmadd(const RealVec & a, const RealVec& b, const RealVec& c ) {
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return c - (a * b);
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fnmadd_pd(a, b, c);
	#else /* KNL */
		return _mm512_fnmadd_pd(a, b, c);
	#endif
	}

	vcp_inline
	static RealVec fmsub(const RealVec & a, const RealVec& b, const RealVec& c) {
	#if VCP_VEC_TYPE == VCP_NOVEC or VCP_VEC_TYPE == VCP_VEC_SSE3 or VCP_VEC_TYPE == VCP_VEC_AVX
		return (a * b) - c;
	#elif VCP_VEC_TYPE == VCP_VEC_AVX2
		return _mm256_fmsub_pd(a, b, c);
	#else /* KNL */
		return _mm512_fmsub_pd(a, b, c);
	#endif
	}

	vcp_inline
	RealVec operator / (const RealVec& rhs) const {
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

	vcp_inline
	static RealVec max (const RealVec& a, const RealVec& b) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return std::max(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_max_pd(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_max_pd(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_max_pd(a, b);
	#endif
	}

	vcp_inline
	static RealVec cos (const RealVec& a) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return std::cos(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_cos_pd(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_cosd_pd(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cosd_pd(a);
	#endif
	}

	vcp_inline
	static RealVec sqrt (const RealVec& rhs) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return std::sqrt(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_sqrt_pd(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_sqrt_pd(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_sqrt_pd(rhs);
	#endif
	}

	vcp_inline
	static RealVec scal_prod(
		const RealVec& a1, const RealVec& a2, const RealVec& a3,
		const RealVec& b1, const RealVec& b2, const RealVec& b3) {
		return fmadd(a1, b1, fmadd(a2, b2, a3 * b3));
	}

	vcp_inline
	static RealVec set1(const double& v) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return v;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_set1_pd(v);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_set1_pd(v);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_set1_pd(v);
	#endif
	}

	vcp_inline
	static RealVec aligned_load(const double * const a) {
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

	vcp_inline
	static RealVec broadcast(const double * const a) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return *a;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_loaddup_pd(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_broadcast_sd(a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_set1_pd(*a);
	#endif
	}

	vcp_inline
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

	vcp_inline
	static RealVec unpack_lo(const RealVec& a, const RealVec& b) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0 / 0.0; // makes no sense
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpacklo_pd(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpacklo_pd(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0 / 0.0); // not used!
	#endif
	}

	vcp_inline
	static RealVec unpack_hi(const RealVec& a, const RealVec& b) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0 / 0.0; // makes no sense
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_unpackhi_pd(a,b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_unpackhi_pd(a, b);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return set1(0.0 / 0.0); // not used!
	#endif
	}

	vcp_inline
	MaskVec<double> operator < (const RealVec & rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d < rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_RealCalcVec_to_MaskVec(_mm_cmplt_pd(_d, rhs));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_RealCalcVec_to_MaskVec(_mm256_cmp_pd(_d, rhs, _CMP_LT_OS));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cmp_pd_mask(_d, rhs, _CMP_LT_OS);
	#endif
	}

	vcp_inline
	MaskVec<double> operator != (const RealVec & rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return _d != rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return cast_RealCalcVec_to_MaskVec(_mm_cmpneq_pd(_d, rhs));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return cast_RealCalcVec_to_MaskVec(_mm256_cmp_pd(_d, rhs, _CMP_NEQ_OS));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cmp_pd_mask(_d, rhs, _CMP_NEQ_UQ);
	#endif
	}

	vcp_inline
	static RealVec apply_mask(const RealVec& d, const MaskVec<double>& m) {
	#if VCP_VEC_TYPE == VCP_NOVEC
		return m ? d : RealVec::zero();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_mov_pd(RealVec::zero(), m, d);
	#else // SSE, AVX, AVX2
		return cast_MaskVec_to_RealCalcVec(cast_RealCalcVec_to_MaskVec(d) and m);
	#endif
	}

	vcp_inline
	static RealVec aligned_load_mask(const double * const a, MaskVec<double> m) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return apply_mask(RealVec::aligned_load(a),m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return apply_mask(RealVec::aligned_load(a),m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_maskload_pd(a, m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_mask_load_pd(RealVec::zero(), m, a);
	#endif
	}

	vcp_inline
	static void horizontal_add_and_store(const RealVec& a, double * const mem_addr) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		// add one double
	    (*mem_addr) += a;

	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		// add two doubles
		const RealVec t1 = _mm_hadd_pd(a, a);
		const RealVec t2 = _mm_add_sd( t1, _mm_load_sd(mem_addr));
		_mm_store_sd( mem_addr, t2);

	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		// add four doubles
		static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
		const RealVec a_t1 = _mm256_permute2f128_pd(a, a, 0x1);
		const RealVec a_t2 = _mm256_hadd_pd(a, a_t1);
		const RealVec a_t3 = _mm256_hadd_pd(a_t2, a_t2);
		_mm256_maskstore_pd(
			mem_addr,
			memoryMask_first,
				a_t3 + RealVec::aligned_load_mask(mem_addr, memoryMask_first)
		);

	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		// add eight doubles
		__m256d low = _mm512_extractf64x4_pd(a, 0);
		__m256d high = _mm512_extractf64x4_pd(a, 1);

		__m256d t1 = _mm256_hadd_pd(low + high, low + high);
		__m128d t2 = _mm256_extractf128_pd(t1,1);
		__m128d t3 = _mm_add_sd(_mm256_castpd256_pd128(t1),t2);

		*mem_addr += _mm_cvtsd_f64(t3);
	#endif
	}

	vcp_inline
	void aligned_load_add_store(double * location) const {
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
	static RealVec fastReciprocal_mask(const RealVec& d, const MaskVec<double>& m) {

		/* reciprocal */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const __m128 denom_ps = _mm256_cvtpd_ps(d);
				const __m128 recip_ps = _mm_rcp_ps(denom_ps);

				const real_vec inv_unmasked = _mm256_cvtps_pd(recip_ps); //12bit
			#elif VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL
				const RealVec inv = _mm512_maskz_rcp28_pd(m, d); //28bit
			#elif VCP_VEC_TYPE == VCP_VEC_AVX512F or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				const RealVec inv = _mm512_maskz_rcp14_pd(m, d); //14bit
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealVec inv_unmasked = set1(1.0) / d;
			#endif

		/* mask and/or N-R-iterations if necessary */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const RealVec inv = apply_mask(inv_unmasked, m); //12bit

				const RealVec inv_24bits = inv * fnmadd(d, inv, set1(2.0)); 				//24bit, 1. N-R-Iteration
				const RealVec inv_48bits = inv_24bits * fnmadd(d, inv_24bits, set1(2.0));	//48bit, 2. N-R-Iteration
				const RealVec inv_prec = inv_48bits * fnmadd(d, inv_48bits, set1(2.0)); 	//96bit, 3. N-R-Iteration
			#elif VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL
				const RealVec inv_prec = inv * fnmadd(d, inv, set1(2.0)); //56bit, 1 N-R-Iteration
			#elif VCP_VEC_TYPE == VCP_VEC_AVX512F or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				const RealVec inv_28bits = inv * fnmadd(d, inv, set1(2.0));					//28bit, 1. N-R-Iteration
				const RealVec inv_prec = inv_28bits * fnmadd(d, inv_28bits, set1(2.0)); 	//56bit, 2. N-R-Iteration
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealVec inv_prec = apply_mask(inv_unmasked, m);
			#endif

		/* check for NaNs */
		#ifndef NDEBUG
			//set nan_found to zero only if NO NaN is found
				#if   VCP_VEC_WIDTH == VCP_VEC_W__64
					int nan_found = std::isnan(inv_prec) ? 1 : 0;
				#elif VCP_VEC_WIDTH == VCP_VEC_W_128
					int nan_found = _mm_movemask_pd(_mm_cmpunord_pd(inv_prec, inv_prec));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_256
					__m128d low = _mm256_castpd256_pd128(inv_prec);
					__m128d high = _mm256_extractf128_pd(inv_prec, 1);
					int nan_found = _mm_movemask_pd(_mm_cmpunord_pd(low, high));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_512
					MaskVec<double> nan_found = _mm512_cmp_pd_mask(inv_prec, inv_prec, _CMP_UNORD_Q);
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
	static RealVec fastReciprocSqrt_mask(const RealVec& d, const MaskVec<double>& m) {

		/* reciprocal sqrt */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const __m128 denom_ps = _mm256_cvtpd_ps(d);
				const __m128 recipSqrt_ps = _mm_rsqrt_ps(denom_ps);

				const real_vec invSqrt_unmasked = _mm256_cvtps_pd(recipSqrt_ps); //12bit
			#elif VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL
				const RealVec invSqrt = _mm512_maskz_rsqrt28_pd(m, d); //28bit
			#elif VCP_VEC_TYPE == VCP_VEC_AVX512F or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				const RealVec invSqrt = _mm512_maskz_rsqrt14_pd(m, d); //14bit
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealVec invSqrt_unmasked = set1(1.0) / sqrt(d);
			#endif

		/* mask and/or N-R-iterations if necessary */
			#if VCP_VEC_WIDTH == VCP_VEC_W_256 and VCP_VEC_TYPE == VCP_VEC_AVX2
				const RealVec invSqrt = apply_mask(invSqrt_unmasked, m); //12bit

				const RealVec d2 = d * set1(0.5);
				const RealVec invSqrt_24bits = invSqrt 		  * fnmadd(d2, invSqrt * invSqrt, set1(1.5)); 				//24bit, 1. N-R-Iteration
				const RealVec invSqrt_48bits = invSqrt_24bits * fnmadd(d2, invSqrt_24bits * invSqrt_24bits, set1(1.5));	//48bit, 2. N-R-Iteration
				const RealVec invSqrt_prec   = invSqrt_48bits * fnmadd(d2, invSqrt_48bits * invSqrt_48bits, set1(1.5)); //96bit, 3. N-R-Iteration
			#elif VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL
				const RealVec invSqrt_prec = invSqrt * fnmadd(d * set1(0.5), invSqrt * invSqrt, set1(1.5)); //56bit, 1 N-R-Iteration
			#elif VCP_VEC_TYPE == VCP_VEC_AVX512F or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
				const RealVec d2 = d * set1(0.5);
				const RealVec invSqrt_28bits = invSqrt        *fnmadd(d2, invSqrt * invSqrt, set1(1.5)); // 28bit, 1 N-R-Iteration
				const RealVec invSqrt_prec = invSqrt_28bits *fnmadd(d2, invSqrt_28bits * invSqrt_28bits, set1(1.5)); // 56bit, 2 N-R-Iteration
			#else /* VCP_VEC_WIDTH == 64/128/256AVX */
				const RealVec invSqrt_prec = apply_mask(invSqrt_unmasked, m);
			#endif

		/* check for NaNs */
		#ifndef NDEBUG
			//set nan_found to zero only if NO NaN is found
				#if   VCP_VEC_WIDTH == VCP_VEC_W__64
					int nan_found = std::isnan(invSqrt_prec) ? 1 : 0;
				#elif VCP_VEC_WIDTH == VCP_VEC_W_128
					int nan_found = _mm_movemask_pd(_mm_cmpunord_pd(invSqrt_prec, invSqrt_prec));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_256
					__m128d low = _mm256_castpd256_pd128(invSqrt_prec);
					__m128d high = _mm256_extractf128_pd(invSqrt_prec, 1);
					int nan_found = _mm_movemask_pd(_mm_cmpunord_pd(low, high));
				#elif VCP_VEC_WIDTH == VCP_VEC_W_512
					MaskVec<double> nan_found = _mm512_cmp_pd_mask(invSqrt_prec, invSqrt_prec, _CMP_UNORD_Q);
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

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVECDOUBLE_H_ */
