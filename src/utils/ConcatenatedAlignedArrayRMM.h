/*
 * ConcatenatedAlignedArrayRMM.h
 *
 *  Created on: Sep 30, 2017
 *      Author: tchipevn
 */

#ifndef UTILS_CONCATENATEDALIGNEDARRAYRMM_H_
#define UTILS_CONCATENATEDALIGNEDARRAYRMM_H_

#include <utils/AlignedArray.h>
#include <particleContainer/adapter/vectorization/SIMD_TYPES.h>
#include <array>

template <typename real_calc_t, typename real_accum_t, typename uid_t>
class ConcatenatedAlignedArrayRMM {
public:
	enum class Quantity_t { RX = 0, RY = 1, RZ = 2, VX = 3, VY = 4, VZ = 5, UID = 6};

	ConcatenatedAlignedArrayRMM(size_t initialSize = 0) : _byteBuffer(0), _numEntriesPerArray(0) {
		mardyn_assert(sizeof(real_calc_t) <= sizeof(uid_t));
		ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::resize(initialSize);
	}

	void prefetchForForce() const;

	real_calc_t* begin_calc(Quantity_t coord);
	real_accum_t* begin_accum(Quantity_t coord);
	uid_t* begin_uid(Quantity_t coord);
	real_calc_t& get_calc(Quantity_t coord, size_t i);
	real_accum_t& get_accum(Quantity_t coord, size_t i);
	uid_t& get_uid(Quantity_t coord, size_t i);

	const real_calc_t* begin_calc(Quantity_t coord) const;
	const real_accum_t* begin_accum(Quantity_t coord) const;
	const uid_t* begin_uid(Quantity_t coord) const;
	const real_calc_t& get_calc(Quantity_t coord, size_t i) const;
	const real_accum_t& get_accum(Quantity_t coord, size_t i) const;
	const uid_t& get_uid(Quantity_t coord, size_t i) const;

	void zero(size_t start_idx = 0);

	/**
	 * \brief Reallocate the array. All content may be lost.
	 */
	void resize(size_t nEntriesPerArray);

	void increaseStorage(size_t oldNumElements, size_t additionalElements);

	void appendValues(std::array<real_calc_t, 3> calcs, std::array<real_accum_t, 3> accums, uid_t uid, size_t oldNumElements);

	size_t get_dynamic_memory() const {
		return _byteBuffer.get_dynamic_memory();
	}

private:
	typedef unsigned char byte_t;
	byte_t* begin(Quantity_t coord);
	const byte_t* begin(Quantity_t coord) const;

	AlignedArray<byte_t, CACHE_LINE_SIZE> _byteBuffer;

	//! how many entries are allocated per real array. The number may be smaller for UID.
	size_t _numEntriesPerArray;
};

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline typename ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::byte_t* ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::begin(Quantity_t coord) {
	byte_t * rx = _byteBuffer;
	byte_t * ry = rx + _numEntriesPerArray * sizeof(real_calc_t);
	byte_t * rz = ry + _numEntriesPerArray * sizeof(real_calc_t);
	byte_t * vx = rz + _numEntriesPerArray * sizeof(real_calc_t);
	byte_t * vy = vx + _numEntriesPerArray * sizeof(real_accum_t);
	byte_t * vz = vy + _numEntriesPerArray * sizeof(real_accum_t);
	byte_t * uid = vz + _numEntriesPerArray * sizeof(real_accum_t);
	byte_t * starts[7] = {rx, ry, rz, vx, vy, vz, uid};
	return _numEntriesPerArray > 0 ? starts[size_t(coord)] : nullptr;
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline const typename ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::byte_t* ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::begin(Quantity_t coord) const {
	const byte_t * rx = _byteBuffer;
	const byte_t * ry = rx + _numEntriesPerArray * sizeof(real_calc_t);
	const byte_t * rz = ry + _numEntriesPerArray * sizeof(real_calc_t);
	const byte_t * vx = rz + _numEntriesPerArray * sizeof(real_calc_t);
	const byte_t * vy = vx + _numEntriesPerArray * sizeof(real_accum_t);
	const byte_t * vz = vy + _numEntriesPerArray * sizeof(real_accum_t);
	const byte_t * uid = vz + _numEntriesPerArray * sizeof(real_accum_t);
	const byte_t * starts[7] = {rx, ry, rz, vx, vy, vz, uid};
	return _numEntriesPerArray > 0 ? starts[size_t(coord)] : nullptr;
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline real_calc_t* ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::begin_calc(Quantity_t coord) {
	mardyn_assert(coord < Quantity_t::VX);
	byte_t * ret = begin(coord);
	return reinterpret_cast<real_calc_t*>(ret);
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline real_accum_t* ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::begin_accum(Quantity_t coord) {
	mardyn_assert(coord >= Quantity_t::VX and coord < Quantity_t::UID);
	byte_t * ret = begin(coord);
	return reinterpret_cast<real_accum_t*>(ret);
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline uid_t* ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::begin_uid(Quantity_t coord) {
	mardyn_assert(coord == Quantity_t::UID);
	byte_t * ret = begin(coord);
	return reinterpret_cast<uid_t*>(ret);
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline real_calc_t& ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::get_calc(Quantity_t coord, size_t i) {
	mardyn_assert(coord < Quantity_t::VX);
	mardyn_assert(i < _numEntriesPerArray);
	byte_t * startByte = begin(coord);
	real_calc_t * startReal = reinterpret_cast<real_calc_t*>(startByte);
	real_calc_t & ret = startReal[i];
	return ret;
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline real_accum_t& ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::get_accum(Quantity_t coord, size_t i) {
	mardyn_assert(coord >= Quantity_t::VX and coord < Quantity_t::UID);
	mardyn_assert(i < _numEntriesPerArray);
	byte_t * startByte = begin(coord);
	real_accum_t * startReal = reinterpret_cast<real_accum_t*>(startByte);
	real_accum_t & ret = startReal[i];
	return ret;
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline uid_t& ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::get_uid(Quantity_t coord, size_t i) {
	mardyn_assert(coord == Quantity_t::UID);
	mardyn_assert(i < _numEntriesPerArray);
	byte_t * startByte = begin(coord);
	uid_t * startUID = reinterpret_cast<uid_t*>(startByte);
	uid_t & ret = startUID[i];
	return ret;
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline const real_calc_t* ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::begin_calc(Quantity_t coord) const {
	mardyn_assert(coord < Quantity_t::VX);
	const byte_t * ret = begin(coord);
	return reinterpret_cast<const real_calc_t*>(ret);
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline const real_accum_t* ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::begin_accum(Quantity_t coord) const {
	mardyn_assert(coord >= Quantity_t::VX and coord < Quantity_t::UID);
	const byte_t * ret = begin(coord);
	return reinterpret_cast<const real_accum_t*>(ret);
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline const uid_t* ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::begin_uid(Quantity_t coord) const {
	mardyn_assert(coord == Quantity_t::UID);
	const byte_t * ret = begin(coord);
	return reinterpret_cast<const uid_t*>(ret);
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline const real_calc_t& ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::get_calc(Quantity_t coord, size_t i) const {
	mardyn_assert(coord < Quantity_t::VX);
	mardyn_assert(i < _numEntriesPerArray);
	const byte_t * startByte = begin(coord);
	const real_calc_t * startReal = reinterpret_cast<const real_calc_t*>(startByte);
	const real_calc_t & ret = startReal[i];
	return ret;
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline const real_accum_t& ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::get_accum(Quantity_t coord, size_t i) const {
	mardyn_assert(coord >= Quantity_t::VX and coord < Quantity_t::UID);
	mardyn_assert(i < _numEntriesPerArray);
	const byte_t * startByte = begin(coord);
	const real_accum_t * startReal = reinterpret_cast<const real_accum_t*>(startByte);
	const real_accum_t & ret = startReal[i];
	return ret;
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline const uid_t& ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::get_uid(Quantity_t coord, size_t i) const {
	mardyn_assert(coord == Quantity_t::UID);
	mardyn_assert(i < _numEntriesPerArray);
	const byte_t * startByte = begin(coord);
	const uid_t * startUID = reinterpret_cast<const uid_t*>(startByte);
	const uid_t & ret = startUID[i];
	return ret;
}

template <typename real_calc_t, typename real_accum_t, typename uid_t>
inline void ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::resize(size_t nEntriesPerArray) {

	if (nEntriesPerArray == 0 and _numEntriesPerArray == 0)
		return;

	_numEntriesPerArray = AlignedArray<real_calc_t, VCP_ALIGNMENT>::_round_up(nEntriesPerArray);
	size_t numBytesForReals = 3 * _numEntriesPerArray * (sizeof(real_calc_t) + sizeof(real_accum_t));
	size_t numBytesForUID = _numEntriesPerArray * sizeof(uid_t);
	size_t totalNumBytes = numBytesForReals + numBytesForUID;

	_byteBuffer.resize(totalNumBytes);

	if (_numEntriesPerArray > nEntriesPerArray) {
		zero(nEntriesPerArray);
	}
}

template<typename real_calc_t, typename real_accum_t, typename uid_t>
inline void ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::zero(size_t start_idx) {
	if (start_idx == 0) {
		_byteBuffer.zero(start_idx);
		return;
	}
	size_t num_to_zero = _numEntriesPerArray - start_idx;
	if (_numEntriesPerArray > 0 and num_to_zero > 0) {
		const int qbegin = static_cast<int>(Quantity_t::RX);
		const int qmid = static_cast<int>(Quantity_t::VX);
		const int qend = static_cast<int>(Quantity_t::UID);
		for (int i = qbegin; i < qmid; ++i) {
			Quantity_t q = static_cast<Quantity_t>(i);
			std::memset(&(get_calc(q,start_idx)), 0, num_to_zero * sizeof(real_calc_t));
		}
		for (int i = qmid; i < qend; ++i) {
			Quantity_t q = static_cast<Quantity_t>(i);
			std::memset(&(get_accum(q,start_idx)), 0, num_to_zero * sizeof(real_accum_t));
		}

		size_t startIndex = begin(Quantity_t::UID) - begin(Quantity_t::RX) + start_idx * sizeof(uid_t);
		_byteBuffer.zero(startIndex);
	}
}

template<typename real_calc_t, typename real_accum_t, typename uid_t>
inline void ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::increaseStorage(size_t oldNumElements, size_t additionalElements) {
	mardyn_assert(oldNumElements <= _numEntriesPerArray);

	size_t newNumElements = oldNumElements + additionalElements;

	if (newNumElements <= _numEntriesPerArray) {
		// no need to resize
		return;
	}

	// do we need to keep contents?
	if (oldNumElements > 0) {
		// yes
		AlignedArray<byte_t, CACHE_LINE_SIZE> backupCopy(_byteBuffer);

		size_t oldNumEntriesPerArray = _numEntriesPerArray;
		resize(newNumElements);

		const int qbegin = static_cast<int>(Quantity_t::RX);
		const int qend = static_cast<int>(Quantity_t::UID);
		for (int i = qbegin; i < qend; ++i) {
			Quantity_t q = static_cast<Quantity_t>(i);
			std::memcpy(begin(q), &(backupCopy[i*oldNumEntriesPerArray*sizeof(real_calc_t)]), oldNumElements * sizeof(real_calc_t));
		}
		std::memcpy(begin(Quantity_t::UID), &(backupCopy[6*oldNumEntriesPerArray*sizeof(real_calc_t)]), oldNumElements * sizeof(uid_t));
	} else {
		// no
		resize(newNumElements);
	}
}

template<typename real_calc_t, typename real_accum_t, typename uid_t>
inline void ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::appendValues(
		std::array<real_calc_t, 3> calcs, std::array<real_accum_t, 3> accums, uid_t uid, size_t oldNumElements) {
	mardyn_assert(oldNumElements <= _numEntriesPerArray);
	if (oldNumElements < _numEntriesPerArray) {
		// no need to resize, baby
	} else {
		increaseStorage(oldNumElements, 1);
	}
	get_calc(Quantity_t::RX, oldNumElements) = calcs[0];
	get_calc(Quantity_t::RY, oldNumElements) = calcs[1];
	get_calc(Quantity_t::RZ, oldNumElements) = calcs[2];
	get_accum(Quantity_t::VX, oldNumElements) = accums[0];
	get_accum(Quantity_t::VY, oldNumElements) = accums[1];
	get_accum(Quantity_t::VZ, oldNumElements) = accums[2];
	get_uid(Quantity_t::UID, oldNumElements) = uid;
}

template<typename real_calc_t, typename real_accum_t, typename uid_t>
inline void ConcatenatedAlignedArrayRMM<real_calc_t, real_accum_t, uid_t>::prefetchForForce() const {
	// prefetch all up to uid begin
	const byte_t * b = begin(Quantity_t::RX);
	const byte_t * e = begin(Quantity_t::UID);
	const size_t stride = CACHE_LINE_SIZE / sizeof(byte_t);

	for (const byte_t * p = b; p < e; p += stride) {
#if defined(__SSE3__)
		_mm_prefetch(reinterpret_cast<const char*>(p), _MM_HINT_T1);
#else
#endif
	}
}

#endif /* UTILS_CONCATENATEDALIGNEDARRAYRMM_H_ */
