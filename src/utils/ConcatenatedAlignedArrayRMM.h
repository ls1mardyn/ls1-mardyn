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

template <typename TypeReal, typename TypeUID>
class ConcatenatedAlignedArrayRMM {
public:
	enum class Quantity_t { RX = 0, RY = 1, RZ = 2, VX = 3, VY = 4, VZ = 5, UID = 6};

	ConcatenatedAlignedArrayRMM(size_t initialSize = 0) : _byteBuffer(0), _numEntriesPerRealArray(0) {
		mardyn_assert(sizeof(TypeReal) <= sizeof(TypeUID));
		ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::resize(initialSize);
	}

	void prefetchForForce() const;

	TypeReal* begin_real(Quantity_t coord);
	TypeUID* begin_uid(Quantity_t coord);
	TypeReal& get_real(Quantity_t coord, size_t i);
	TypeUID& get_uid(Quantity_t coord, size_t i);

	const TypeReal*  begin_real(Quantity_t coord) const;
	const TypeUID* begin_uid(Quantity_t coord) const;
	const TypeReal& get_real(Quantity_t coord, size_t i) const;
	const TypeUID& get_uid(Quantity_t coord, size_t i) const;

	void zero(size_t start_idx = 0);

	/**
	 * \brief Reallocate the array. All content may be lost.
	 */
	void resize(size_t nEntriesPerArray);

	void increaseStorage(size_t oldNumElements, size_t additionalElements);

	void appendValues(std::array<TypeReal,6> reals, TypeUID uid, size_t oldNumElements);

	size_t get_dynamic_memory() const {
		return _byteBuffer.get_dynamic_memory();
	}

private:
	typedef unsigned char byte_t;
	byte_t* begin(Quantity_t coord);
	const byte_t* begin(Quantity_t coord) const;

	AlignedArray<byte_t, CACHE_LINE_SIZE> _byteBuffer;

	//! how many entries are allocated per real array. The number may be smaller for UID.
	size_t _numEntriesPerRealArray;
};

template <typename TypeReal, typename TypeUID>
inline typename ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::byte_t* ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::begin(Quantity_t coord) {
	byte_t * data = _byteBuffer;
	return _numEntriesPerRealArray > 0 ?
			data + (size_t(coord) * _numEntriesPerRealArray * sizeof(TypeReal)) : nullptr;
}

template <typename TypeReal, typename TypeUID>
inline const typename ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::byte_t* ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::begin(Quantity_t coord) const {
	const byte_t * data = _byteBuffer;
	return _numEntriesPerRealArray > 0 ?
			data + (size_t(coord) * _numEntriesPerRealArray * sizeof(TypeReal)) : nullptr;
}

template <typename TypeReal, typename TypeUID>
inline TypeReal* ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::begin_real(Quantity_t coord) {
	mardyn_assert(coord != Quantity_t::UID);
	byte_t * ret = begin(coord);
	return reinterpret_cast<TypeReal*>(ret);
}

template <typename TypeReal, typename TypeUID>
inline TypeUID* ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::begin_uid(Quantity_t coord) {
	mardyn_assert(coord == Quantity_t::UID);
	byte_t * ret = begin(coord);
	return reinterpret_cast<TypeUID*>(ret);
}

template <typename TypeReal, typename TypeUID>
inline TypeReal& ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::get_real(Quantity_t coord, size_t i) {
	mardyn_assert(coord != Quantity_t::UID);
	mardyn_assert(i < _numEntriesPerRealArray);
	byte_t * startByte = begin(coord);
	TypeReal * startReal = reinterpret_cast<TypeReal*>(startByte);
	TypeReal & ret = startReal[i];
	return ret;
}

template <typename TypeReal, typename TypeUID>
inline TypeUID& ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::get_uid(Quantity_t coord, size_t i) {
	mardyn_assert(coord == Quantity_t::UID);
	mardyn_assert(i < _numEntriesPerRealArray);
	byte_t * startByte = begin(coord);
	TypeUID * startUID = reinterpret_cast<TypeUID*>(startByte);
	TypeUID & ret = startUID[i];
	return ret;
}

template <typename TypeReal, typename TypeUID>
inline const TypeReal* ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::begin_real(Quantity_t coord) const {
	mardyn_assert(coord != Quantity_t::UID);
	const byte_t * ret = begin(coord);
	return reinterpret_cast<const TypeReal*>(ret);
}

template <typename TypeReal, typename TypeUID>
inline const TypeUID* ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::begin_uid(Quantity_t coord) const {
	mardyn_assert(coord == Quantity_t::UID);
	const byte_t * ret = begin(coord);
	return reinterpret_cast<const TypeUID*>(ret);
}

template <typename TypeReal, typename TypeUID>
inline const TypeReal& ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::get_real(Quantity_t coord, size_t i) const {
	mardyn_assert(coord != Quantity_t::UID);
	mardyn_assert(i < _numEntriesPerRealArray);
	const byte_t * startByte = begin(coord);
	const TypeReal * startReal = reinterpret_cast<const TypeReal*>(startByte);
	const TypeReal & ret = startReal[i];
	return ret;
}

template <typename TypeReal, typename TypeUID>
inline const TypeUID& ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::get_uid(Quantity_t coord, size_t i) const {
	mardyn_assert(coord == Quantity_t::UID);
	mardyn_assert(i < _numEntriesPerRealArray);
	const byte_t * startByte = begin(coord);
	const TypeUID * startUID = reinterpret_cast<const TypeUID*>(startByte);
	const TypeUID & ret = startUID[i];
	return ret;
}

template <typename TypeReal, typename TypeUID>
inline void ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::resize(size_t nEntriesPerArray) {

	if (nEntriesPerArray == 0 and _numEntriesPerRealArray == 0)
		return;

	_numEntriesPerRealArray = AlignedArray<TypeReal, VCP_ALIGNMENT>::_round_up(nEntriesPerArray);
	size_t numBytesForReals = 6 * _numEntriesPerRealArray * sizeof(TypeReal);
	size_t numBytesForUID = AlignedArray<TypeUID, VCP_ALIGNMENT>::_round_up(nEntriesPerArray) * sizeof(TypeUID); //TODO: VCP_ALIGNMENT should not be needed here
	size_t totalNumBytes = numBytesForReals + numBytesForUID;

	_byteBuffer.resize(totalNumBytes);

	if (_numEntriesPerRealArray > nEntriesPerArray) {
		zero(nEntriesPerArray);
	}
}

template<typename TypeReal, typename TypeUID>
inline void ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::zero(size_t start_idx) {
	if (start_idx == 0) {
		_byteBuffer.zero(start_idx);
		return;
	}
	size_t num_to_zero = AlignedArray<TypeReal, VCP_ALIGNMENT>::_round_up(start_idx) - start_idx;
	if (_numEntriesPerRealArray > 0 and num_to_zero > 0) {
		const int qbegin = static_cast<int>(Quantity_t::RX);
		const int qend = static_cast<int>(Quantity_t::UID);
		for (int i = qbegin; i < qend; ++i) {
			Quantity_t q = static_cast<Quantity_t>(i);
			std::memset(&(get_real(q,start_idx)), 0, num_to_zero * sizeof(TypeReal));
		}

		size_t startIndex = begin(Quantity_t::UID) - begin(Quantity_t::RX) + start_idx * sizeof(TypeUID);
		_byteBuffer.zero(startIndex);
	}
}

template<typename TypeReal, typename TypeUID>
inline void ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::increaseStorage(size_t oldNumElements, size_t additionalElements) {
	mardyn_assert(oldNumElements <= _numEntriesPerRealArray);

	size_t newNumElements = oldNumElements + additionalElements;

	if (newNumElements <= _numEntriesPerRealArray) {
		// no need to resize
		return;
	}

	// do we need to keep contents?
	if (oldNumElements > 0) {
		// yes
		AlignedArray<byte_t, CACHE_LINE_SIZE> backupCopy(_byteBuffer);

		size_t oldNumEntriesPerArray = _numEntriesPerRealArray;
		resize(newNumElements);

		const int qbegin = static_cast<int>(Quantity_t::RX);
		const int qend = static_cast<int>(Quantity_t::UID);
		for (int i = qbegin; i < qend; ++i) {
			Quantity_t q = static_cast<Quantity_t>(i);
			std::memcpy(begin(q), &(backupCopy[i*oldNumEntriesPerArray*sizeof(TypeReal)]), oldNumElements * sizeof(TypeReal));
		}
		std::memcpy(begin(Quantity_t::UID), &(backupCopy[6*oldNumEntriesPerArray*sizeof(TypeReal)]), oldNumElements * sizeof(TypeUID));
	} else {
		// no
		resize(newNumElements);
	}
}

template<typename TypeReal, typename TypeUID>
inline void ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::appendValues(
		std::array<TypeReal, 6> reals, TypeUID uid, size_t oldNumElements) {
	mardyn_assert(oldNumElements <= _numEntriesPerRealArray);
	if (oldNumElements < _numEntriesPerRealArray) {
		// no need to resize, baby
	} else {
		increaseStorage(oldNumElements, 1);
	}
	const int qbegin = static_cast<int>(Quantity_t::RX);
	const int qend = static_cast<int>(Quantity_t::UID);
	for (int i = qbegin; i < qend; ++i) {
		Quantity_t q = static_cast<Quantity_t>(i);
		get_real(q,oldNumElements) = reals[i];
	}
	get_uid(Quantity_t::UID, oldNumElements) = uid;
}

template<typename TypeReal, typename TypeUID>
inline void ConcatenatedAlignedArrayRMM<TypeReal, TypeUID>::prefetchForForce() const {
	// prefetch all up to uid begin
	const byte_t * b = begin(Quantity_t::RX);
	const byte_t * e = begin(Quantity_t::UID);
	const size_t stride = CACHE_LINE_SIZE / sizeof(byte_t);

	for (const byte_t * p = b; p < e; p += stride) {
#if defined(__SSE3__)
		_mm_prefetch(reinterpret_cast<const char*>(p), _MM_HINT_T1);
#elif defined(__MIC__)
		_mm_prefetch(reinterpret_cast<const char*>(p), 2);
#else
#endif
	}
}

#endif /* UTILS_CONCATENATEDALIGNEDARRAYRMM_H_ */
