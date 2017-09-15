/*
 * AlignedArrayTriplet.h
 *
 *  Created on: Mar 23, 2016
 *      Author: tchipevn
 */

#ifndef ALIGNEDARRAYTRIPLET_H
#define ALIGNEDARRAYTRIPLET_H

#include "AlignedArray.h"
#include "utils/mardyn_assert.h"

template <class T>
class AlignedArrayTriplet : public AlignedArray<T> {
public:
	AlignedArrayTriplet(size_t initialSize = 0) : AlignedArray<T>(0), _numEntriesPerArray(0) {
		AlignedArrayTriplet<T>::resize(initialSize);
	}

	void prefetch(int hint = 1, int n = -1) const {
		mardyn_assert(n >= -2);

		int endPrefetch;
		const int stride = this->_round_up(1);

		switch(n) {
		case -1:
			// prefetch all up to capacity()
			AlignedArray<T>::prefetch(hint, n);
			return;
			break;
		case -2:
			// prefetch all up to size()
			endPrefetch = _numEntriesPerArray;
			break;
		default:
			// prefetch only first n elements
			endPrefetch = n;
		}

		for (int i = 0; i < endPrefetch; i+= stride) {
#if defined(__SSE3__)
			_mm_prefetch((const char*)&(x(i)), _MM_HINT_T1);
			_mm_prefetch((const char*)&(y(i)), _MM_HINT_T1);
			_mm_prefetch((const char*)&(z(i)), _MM_HINT_T1);
#elif defined(__MIC__)
			_mm_prefetch((const char*)&(x(i)), 2);
			_mm_prefetch((const char*)&(y(i)), 2);
			_mm_prefetch((const char*)&(z(i)), 2);
#else
#endif
		}
	}

	T* xBegin() { return _numEntriesPerArray > 0 ? this->_vec.data() + (0 * _numEntriesPerArray) : nullptr; } // TODO: remove ternary operator?
	T* yBegin() { return _numEntriesPerArray > 0 ? this->_vec.data() + (1 * _numEntriesPerArray) : nullptr; }
	T* zBegin() { return _numEntriesPerArray > 0 ? this->_vec.data() + (2 * _numEntriesPerArray) : nullptr; }
	const T* xBegin() const { return _numEntriesPerArray > 0 ? this->_vec.data() + (0 * _numEntriesPerArray) : nullptr; }
	const T* yBegin() const { return _numEntriesPerArray > 0 ? this->_vec.data() + (1 * _numEntriesPerArray) : nullptr; }
	const T* zBegin() const { return _numEntriesPerArray > 0 ? this->_vec.data() + (2 * _numEntriesPerArray) : nullptr; }

	T& x(size_t i) { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 0 * _numEntriesPerArray]; }
	T& y(size_t i) { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 1 * _numEntriesPerArray]; }
	T& z(size_t i) { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 2 * _numEntriesPerArray]; }
	const T& x(size_t i) const { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 0 * _numEntriesPerArray]; }
	const T& y(size_t i) const { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 1 * _numEntriesPerArray]; }
	const T& z(size_t i) const { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 2 * _numEntriesPerArray]; }

	size_t dimensionToOffset(int i) const {
		mardyn_assert(i >= 0 and i < 3);
		return i * _numEntriesPerArray;
	}

	T& linearCrossAccess(size_t i) { mardyn_assert(i < 3*_numEntriesPerArray); return this->_vec[i];}

	size_t resize_zero_shrink(size_t exact_size, bool zero_rest_of_CL = false, bool allow_shrink = false) {
		size_t size_rounded_up = this->_round_up(exact_size);
		size_t size_rounded_up_x3 = size_rounded_up * 3;
		_numEntriesPerArray = size_rounded_up;

		bool need_resize = size_rounded_up_x3 > this->_vec.size() or (allow_shrink and size_rounded_up_x3 < this->_vec.size());

		if (need_resize) {
			AlignedArray<T>::resize(size_rounded_up_x3);
			// resize zero-s all
		} else {
			// we didn't resize, but we might still need to zero the rest of the Cache Line
			if (zero_rest_of_CL) {
				zero(exact_size);
			}
		}

		return _numEntriesPerArray;
	}

	void zero(size_t start_idx) {
		size_t num_to_zero = this->_round_up(start_idx) - start_idx;
		if (_numEntriesPerArray > 0 and num_to_zero > 0) {
			std::memset(&(x(start_idx)), 0, num_to_zero * sizeof(T));
			std::memset(&(y(start_idx)), 0, num_to_zero * sizeof(T));
			std::memset(&(z(start_idx)), 0, num_to_zero * sizeof(T));
		}
	}

	/**
	 * \brief Reallocate the array. All content may be lost.
	 */
	void resize(size_t nEntriesPerArray) {
		if (nEntriesPerArray == 0 and _numEntriesPerArray == 0)
			return;

		// make sure that this is divisible by eight, otherwise we break alignment
		// get the rounded-up number of entries
		_numEntriesPerArray = this->_round_up(nEntriesPerArray);
		AlignedArray<T>::resize(3 * _numEntriesPerArray);

		if (_numEntriesPerArray > nEntriesPerArray) {
			// set elements from nEntriesPerArray to _numEntriesPerArray to zero:
			size_t elements = _numEntriesPerArray - nEntriesPerArray;
			std::memset(&(x(nEntriesPerArray)), 0, elements * sizeof(T));
			std::memset(&(y(nEntriesPerArray)), 0, elements * sizeof(T));
			std::memset(&(z(nEntriesPerArray)), 0, elements * sizeof(T));
		}
	}

	void increaseStorage(size_t oldNumElements, size_t additionalElements) {
		mardyn_assert(oldNumElements <= _numEntriesPerArray);

		size_t newNumElements = oldNumElements + additionalElements;

		if (newNumElements <= _numEntriesPerArray) {
			// no need to resize
			return;
		}

		// do we need to keep contents?
		if (oldNumElements > 0) {
			// yes
			AlignedArray<T> backupCopy(*this);

			size_t oldNumEntriesPerArray = _numEntriesPerArray;
			resize_zero_shrink(newNumElements);

			std::memcpy(xBegin(), &(backupCopy[0*oldNumEntriesPerArray]), oldNumElements * sizeof(T));
			std::memcpy(yBegin(), &(backupCopy[1*oldNumEntriesPerArray]), oldNumElements * sizeof(T));
			std::memcpy(zBegin(), &(backupCopy[2*oldNumEntriesPerArray]), oldNumElements * sizeof(T));
		} else {
			// no
			resize_zero_shrink(newNumElements);
		}
	}

	void appendValueTriplet(T v0, T v1, T v2, size_t oldNumElements) {
		mardyn_assert(oldNumElements <= _numEntriesPerArray);
		if (oldNumElements < _numEntriesPerArray) {
			// no need to resize, baby
		} else {
			increaseStorage(oldNumElements, 1);
		}
		x(oldNumElements) = v0;
		y(oldNumElements) = v1;
		z(oldNumElements) = v2;
	}

private:
	//! at any point, _data contains precisely three times this storage
	size_t _numEntriesPerArray;
};

#endif /* ALIGNEDARRAYTRIPLET_H */

