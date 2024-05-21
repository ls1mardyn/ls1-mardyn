/*
 * MaskVecDouble.h
 *
 *  Created on: 6 Feb 2018
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_MASKVECDOUBLE_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_MASKVECDOUBLE_H_

#include "MaskVec.h"

// keep this file and MaskVecFloat as close as possible, so that they can be examined via diff!

namespace vcp {

template<>
class MaskVec<double> {
private:
	// own typedefs necessary
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		typedef uint8_t mask_vec;
		typedef uint8_t mask_single;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		typedef __m128i mask_vec;
		typedef uint64_t mask_single;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		typedef __m256i mask_vec;
		typedef uint64_t mask_single;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512

		// these can't be put in std::conditional, because it would just be too nice. warnings.
		typedef __mmask8 mask_vec;
		typedef __mmask8 mask_single;

		#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
			typedef __m512i lookupOrMask_vec;
			typedef countertype32 lookupOrMask_single;
		#endif
	#endif

	mask_vec _m;

public:
	vcp_inline
	MaskVec() {}

	vcp_inline
	operator mask_vec() const {
		return _m;
	}

	vcp_inline
	MaskVec(const mask_vec & m) {
		_m = m;
	}

	vcp_inline
	static MaskVec zero() {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return 0;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			return _mm_setzero_si128();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			return _mm256_setzero_si256();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return 0x00;
	#endif
	}

	vcp_inline
	static MaskVec ones() {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return ~0;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			return _mm_set_epi32(~0, ~0, ~0, ~0);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return 0xFF;
	#endif
	}

	vcp_inline
	MaskVec operator and (const MaskVec& rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return _m & rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			return _mm_and_si128(_m, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			return _mm256_and_si256(_m, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return _m & rhs;
	#endif
	}

	vcp_inline
	MaskVec operator or (const MaskVec& rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return _m | rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			return _mm_or_si128(_m, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			return _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(_m), _mm256_castsi256_pd(rhs)));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return _m | rhs;
	#endif
	}

	vcp_inline
	MaskVec operator xor (const MaskVec & rhs) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return _m xor rhs;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			return _mm_xor_si128(_m, rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			return _mm256_castpd_si256(_mm256_xor_pd(_mm256_castsi256_pd(_m), _mm256_castsi256_pd(rhs)));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return _m ^ rhs;
	#endif
	}

	vcp_inline
	static MaskVec aligned_load(const mask_single * const a) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return *a;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			return _mm_load_si128((const __m128i*)a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			return _mm256_load_si256((const __m256i*)a);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return *a; // is this used ?
	#endif
	}


#if VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
	vcp_inline
	static lookupOrMask_vec aligned_load(const lookupOrMask_single * const a) {
		return _mm512_load_epi64(a);
	}
#endif

	vcp_inline
	void aligned_store(mask_single * location) const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			*location = _m;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			_mm_store_si128((__m128i*)location, _m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			_mm256_store_si256((__m256i*)location, _m);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			*location = _m; // is this used ?
	#endif
	}

	vcp_inline
	int movemask() const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return _m != MaskVec::zero();
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
			return _mm_movemask_pd(_mm_castsi128_pd(_m));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
			return _mm256_movemask_pd(_mm256_castsi256_pd(_m));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return _m != MaskVec::zero();
	#endif
	}

	vcp_inline
	int countUnmasked() const {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
			return _m;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128 or VCP_VEC_WIDTH == VCP_VEC_W_256
			return __builtin_popcount(movemask());
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
			return __builtin_popcount(_m);
	#endif
	}
};

} /* namespace vcp */

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_MASKVECDOUBLE_H_ */
