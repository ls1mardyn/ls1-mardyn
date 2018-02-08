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
	RealAccumVecBackend() {}

	static RealAccumVecBackend convertCalcToAccum(const RealCalcVec & rcv) {
		RealVec<double> first = convert_low(rcv);
		RealVec<double> second = convert_high(rcv);
		return RealAccumVecBackend(first, second);
	}

	RealAccumVecBackend(const RealAccumVecBackend& rhs) {
		_first = rhs._first;
		_second = rhs._second;
	}

	RealAccumVecBackend(const RealVec<double>& first, const RealVec<double>& second) {
		_first = first;
		_second = second;
	}

	RealAccumVecBackend& operator=(const RealAccumVecBackend& rhs) {
		_first = rhs._first;
		_second = rhs._second;
		return *this;
	}

	static RealAccumVecBackend zero() {
		RealAccumVecBackend result;
		result._first = RealVec<double>::zero();
		result._second = RealVec<double>::zero();
		return result;
	}

	RealAccumVecBackend operator+(const RealAccumVecBackend& rhs) const {
		RealAccumVecBackend result;
		result._first = _first + rhs._first;
		result._second = _second + rhs._second;
		return result;
	}

	RealAccumVecBackend operator*(const RealAccumVecBackend& rhs) const {
		RealAccumVecBackend result;
		result._first = _first * rhs._first;
		result._second = _second * rhs._second;
		return result;
	}

	RealAccumVecBackend operator-(const RealAccumVecBackend& rhs) const {
		RealAccumVecBackend result;
		result._first = _first - rhs._first;
		result._second = _second - rhs._second;
		return result;
	}

	static RealAccumVecBackend fmadd(const RealAccumVecBackend & a, const RealAccumVecBackend& b, const RealAccumVecBackend& c ) {
		RealAccumVecBackend result;
		result._first = RealVec<double>::fmadd(a._first, b._first, c._first);
		result._second = RealVec<double>::fmadd(a._second, b._second, c._second);
		return result;
	}

	static RealAccumVecBackend fnmadd(const RealAccumVecBackend & a, const RealAccumVecBackend& b, const RealAccumVecBackend& c ) {
		RealAccumVecBackend result;
		result._first = RealVec<double>::fnmadd(a._first, b._first, c._first);
		result._second = RealVec<double>::fnmadd(a._second, b._second, c._second);
		return result;
	}

	void aligned_store(double * location) const {
		const size_t offset = sizeof(RealVec<double>) / sizeof(double);
		_first.aligned_store(location);
		_second.aligned_store(location + offset);
	}

	static RealAccumVecBackend aligned_load(const double * const a) {
		const size_t offset = sizeof(RealVec<double>) / sizeof(double);
		RealVec<double> first = RealVec<double>::aligned_load(a);
		RealVec<double> second = RealVec<double>::aligned_load(a + offset);
		return RealAccumVecBackend(first, second);
	}

	static RealAccumVecBackend set1(const double& v) {
		RealVec<double> first = RealVec<double>::set1(v);
		RealVec<double> second = RealVec<double>::set1(v);
		return RealAccumVecBackend(first, second);
	}

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
};

#elif VCP_VEC_WIDTH == VCP_VEC_W__64

class RealAccumVecBackend : public RealVec<double> {
public:
	RealAccumVecBackend() {}
	RealAccumVecBackend(const RealVec<double>& rcv) : RealVec<double>() {
		this->_d = rcv;
	}
	static RealAccumVecBackend convertCalcToAccum(const RealCalcVec & rcv) {
		RealAccumVecBackend result;
		result._d = rcv;
		return result;
	}
};

#endif /* VCP_VEC_WIDTH */

} /* namespace vcp */

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALACCUMVECBACKEND_H_ */
