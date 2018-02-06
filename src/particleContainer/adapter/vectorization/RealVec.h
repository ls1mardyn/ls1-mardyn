/*
 * RealVec.h
 *
 *  Created on: 10 Nov 2016
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVEC_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVEC_H_

#include "SIMD_TYPES.h"
#include "MaskVec.h"
#include "utils/mardyn_assert.h"
#include <fstream>
#include <cmath>

#include "utils/Logger.h"

namespace vcp {

template<typename FloatOrDouble>
class RealVec {
private:

public:
	RealVec() {}

	static RealVec cast_MaskVec_to_RealCalcVec(const MaskVec<FloatOrDouble>& m) {
		return RealVec();
	}

	static MaskVec<FloatOrDouble> cast_RealCalcVec_to_MaskVec(const RealVec& d) {
		return RealVec();
	}

	static RealVec zero() {
		return RealVec();
	}

	static RealVec ones() {
		return RealVec();
	}

	RealVec operator + (const RealVec& rhs) const {
		return RealVec();
	}

	RealVec operator - (const RealVec& rhs) const {
		return RealVec();
	}

	RealVec operator * (const RealVec& rhs) const {
		return RealVec();
	}

	static RealVec fmadd(const RealVec & a, const RealVec& b, const RealVec& c ) {
		return RealVec();
	}

	static RealVec fnmadd(const RealVec & a, const RealVec& b, const RealVec& c ) {
		return RealVec();
	}

	static RealVec fmsub(const RealVec & a, const RealVec& b, const RealVec& c) {
		return RealVec();
	}

	RealVec operator / (const RealVec& rhs) const {
		return RealVec();
	}

	static RealVec sqrt (const RealVec& rhs) {
		return RealVec();
	}

	static RealVec scal_prod(
		const RealVec& a1, const RealVec& a2, const RealVec& a3,
		const RealVec& b1, const RealVec& b2, const RealVec& b3) {
		return RealVec();
	}

	static RealVec set1(const vcp_real_calc& v) {
		return RealVec();
	}

	static RealVec aligned_load(const vcp_real_calc * const a) {
		return RealVec();
	}

	static RealVec broadcast(const vcp_real_calc * const a) {
		return RealVec();
	}

	void aligned_store(vcp_real_calc * location) const {

	}

	static RealVec unpack_lo(const RealVec& a, const RealVec& b) {
		return RealVec();
	}

	static RealVec unpack_hi(const RealVec& a, const RealVec& b) {
		return RealVec();
	}

	MaskVec<FloatOrDouble> operator < (const RealVec & rhs) const {
		return MaskVec<FloatOrDouble>();
	}

	MaskVec<FloatOrDouble> operator != (const RealVec & rhs) const {
		return MaskVec<FloatOrDouble>();
	}

	static RealVec apply_mask(const RealVec& d, const MaskVec<FloatOrDouble>& m) {
		return RealVec();
	}

	static RealVec aligned_load_mask(const vcp_real_calc * const a, MaskVec<FloatOrDouble> m) {
		return RealVec();
	}

	vcp_inline
	static void horizontal_add_and_store(const RealVec& a, vcp_real_calc * const mem_addr) {

	}

	void aligned_load_add_store(vcp_real_calc * location) const {
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
	// KNC - small speed-up, because _mm512_div actually already uses Newton-Raphson, but doesn't assume that conversions double -> float are safe?
	// KNL - an educated guess would assume that AVX512ER is there for a reason :)
	static RealVec fastReciprocal_mask(const RealVec& d, const MaskVec<FloatOrDouble>& m) {
		return RealVec();
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
	static RealVec fastReciprocSqrt_mask(const RealVec& d, const MaskVec<FloatOrDouble>& m) {
		return RealVec();
	} //fastReciprocSqrt_mask(..)

}; /* class RealVec */

} /* namespace vcp */

// SPECIALIZATIONS
#include "RealVecFloat.h"
#include "RealVecDouble.h"



#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVEC_H_ */
