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

// this is a somewhat strange construction:
// we use BOTH inheritance
// AND delegation
// from the SAME class
// which needs to work with RealCalcVec, which may be a different class!

template<typename FloatOrDouble>
class RealAccumVecBackend : public RealVec<FloatOrDouble> {

#if VCP_PREC == VCP_SPDP and VCP_VEC_WIDTH != VCP_VEC_W__64
private:
	RealVec<FloatOrDouble> _second;
#endif

public:

	RealAccumVecBackend() : RealVec<FloatOrDouble>() {}

#if not(VCP_PREC == VCP_SPDP and VCP_VEC_WIDTH != VCP_VEC_W__64)
	RealAccumVecBackend(const RealCalcVec & rcv) : RealVec<FloatOrDouble>() {
		this->_d = rcv;
	}

	RealAccumVecBackend& operator=(const RealVec<FloatOrDouble>& rhs) {
		this->_d = rhs;
		return *this;
	}
#endif

#if VCP_PREC == VCP_SPDP and VCP_VEC_WIDTH == VCP_VEC_W__64
	RealAccumVecBackend(const RealVec<FloatOrDouble>& rcv) : RealVec<FloatOrDouble>() {
		this->_d = rcv;
	}
#endif



#if VCP_PREC == VCP_SPDP and VCP_VEC_WIDTH != VCP_VEC_W__64

public:

	RealAccumVecBackend(const RealCalcVec & rcv) : RealVec<FloatOrDouble>() {
		this->_d = convert_low(rcv);
		_second = convert_high(rcv);
	}

	RealAccumVecBackend(const RealAccumVecBackend& rhs) : RealVec<FloatOrDouble>() {
		this->_d = rhs._d;
		_second = rhs._second;
	}

	RealAccumVecBackend(const RealVec<FloatOrDouble>& first, const RealVec<FloatOrDouble>& second) : RealVec<FloatOrDouble>() {
		this->_d = first;
		_second = second;
	}

	RealAccumVecBackend& operator=(const RealAccumVecBackend& rhs) {
		this->_d = rhs._d;
		_second = rhs._second;
		return *this;
	}

	static RealAccumVecBackend zero() {
		RealAccumVecBackend result;
		result._d = RealVec<FloatOrDouble>::zero();
		result._second = RealVec<FloatOrDouble>::zero();
		return result;
	}

	RealAccumVecBackend operator+(const RealAccumVecBackend& rhs) const {
		RealAccumVecBackend result;
		result._d = this->_d + rhs._d;
		result._second = this->_second + rhs._second;
		return result;
	}

	RealAccumVecBackend operator-(const RealAccumVecBackend& rhs) const {
		RealAccumVecBackend result;
		result._d = this->_d - rhs._d;
		result._second = this->_second - rhs._second;
		return result;
	}

	static RealAccumVecBackend fmadd(const RealAccumVecBackend & a, const RealAccumVecBackend& b, const RealAccumVecBackend& c ) {
		RealAccumVecBackend result;
		result._d = fmadd(a._d, b._d, c._d);
		result._second = fmadd(a._second, b._second, c._second);
		return result;
	}

	static RealAccumVecBackend fnmadd(const RealAccumVecBackend & a, const RealAccumVecBackend& b, const RealAccumVecBackend& c ) {
		RealAccumVecBackend result;
		result._d = fnmadd(a._d, b._d, c._d);
		result._second = fnmadd(a._second, b._second, c._second);
		return result;
	}

	void aligned_store(FloatOrDouble * location) const {
		RealVec<FloatOrDouble>::aligned_store(location);
		const size_t offset = sizeof(RealVec<FloatOrDouble>) / sizeof(FloatOrDouble);
		_second.aligned_store(location + offset);
	}

	static RealAccumVecBackend aligned_load(const FloatOrDouble * const a) {
		RealVec<FloatOrDouble> first = RealVec<FloatOrDouble>::aligned_load(a);
		const size_t offset = sizeof(RealVec<FloatOrDouble>) / sizeof(FloatOrDouble);
		RealVec<FloatOrDouble> second = RealVec<FloatOrDouble>::aligned_load(a + offset);
		return RealAccumVecBackend(first, second);
	}

	static RealVec<FloatOrDouble> convert_low(const RealCalcVec& rhs) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		// the NOVEC case is not all that relevant and hard to handle here
		// use both functions
		return static_cast<vcp_real_accum>(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_cvtps_pd(rhs);
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_cvtps_pd(_mm256_extractf128_ps(rhs, 0));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(rhs), 0)));
	#endif

	}

	static RealVec<FloatOrDouble> convert_high(const RealCalcVec& rhs) {
	#if   VCP_VEC_WIDTH == VCP_VEC_W__64
		return 0.0;
	#elif VCP_VEC_WIDTH == VCP_VEC_W_128
		return _mm_cvtps_pd(_mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(rhs), 8)));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_256
		return _mm256_cvtps_pd(_mm256_extractf128_ps(rhs, 1));
	#elif VCP_VEC_WIDTH == VCP_VEC_W_512
		return _mm512_cvtps_pd(_mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(rhs), 1)));
	#endif
	}


#endif /* VCP_SPDP */
};

} /* namespace vcp */

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALACCUMVECBACKEND_H_ */
