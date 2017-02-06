/*
 * AlignedArrayTriplet.h
 *
 *  Created on: Mar 23, 2016
 *      Author: tchipevn
 */

#ifndef ALIGNEDARRAYTRIPLET_H
#define ALIGNEDARRAYTRIPLET_H

#include "AlignedArray.h"
#include <cassert>

template <class T>
class AlignedArrayTriplet : public AlignedArray<T> {
public:
	AlignedArrayTriplet(size_t initialSize = 0) : AlignedArray<T>(0), _numEntriesPerArray(0) {
		AlignedArrayTriplet<T>::resize(initialSize);
	}

	T* xBegin() { return _numEntriesPerArray > 0 ? this->_p + (0 * _numEntriesPerArray) : nullptr; }
	T* yBegin() { return _numEntriesPerArray > 0 ? this->_p + (1 * _numEntriesPerArray) : nullptr; }
	T* zBegin() { return _numEntriesPerArray > 0 ? this->_p + (2 * _numEntriesPerArray) : nullptr; }
	T* xBegin() const { return _numEntriesPerArray > 0 ? this->_p + (0 * _numEntriesPerArray) : nullptr; }
	T* yBegin() const { return _numEntriesPerArray > 0 ? this->_p + (1 * _numEntriesPerArray) : nullptr; }
	T* zBegin() const { return _numEntriesPerArray > 0 ? this->_p + (2 * _numEntriesPerArray) : nullptr; }

	T& x(size_t i) { assert(i < _numEntriesPerArray); return this->_p[i + 0 * _numEntriesPerArray]; }
	T& y(size_t i) { assert(i < _numEntriesPerArray); return this->_p[i + 1 * _numEntriesPerArray]; }
	T& z(size_t i) { assert(i < _numEntriesPerArray); return this->_p[i + 2 * _numEntriesPerArray]; }
	T& x(size_t i) const { assert(i < _numEntriesPerArray); return this->_p[i + 0 * _numEntriesPerArray]; }
	T& y(size_t i) const { assert(i < _numEntriesPerArray); return this->_p[i + 1 * _numEntriesPerArray]; }
	T& z(size_t i) const { assert(i < _numEntriesPerArray); return this->_p[i + 2 * _numEntriesPerArray]; }

	size_t dimensionToOffset(int i) const {
		assert(i >= 0);
		assert(i < 3);
		static size_t rets[3] = {0 * _numEntriesPerArray, 1 * _numEntriesPerArray, 2 * _numEntriesPerArray};
		return rets[i];
	}

	T& linearCrossAccess(size_t i) { assert(i < 3*_numEntriesPerArray); return this->_p[i];}

	size_t resize_zero_shrink(size_t exact_size, bool zero_rest_of_CL = false, bool allow_shrink = false) {
		size_t size_rounded_up = this->_round_up(exact_size);
		size_t size_rounded_up_x3 = size_rounded_up * 3;
		_numEntriesPerArray = size_rounded_up;

		bool need_resize = size_rounded_up_x3 > this->_capacity or (allow_shrink and size_rounded_up_x3 < this->_capacity);

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

	void appendValueTriplet(T v0, T v1, T v2, size_t oldNumElements) {
		assert(oldNumElements <= _numEntriesPerArray);
		if (oldNumElements < _numEntriesPerArray) {
			// no need to resize
		} else {
			// shit, we need to resize, but also keep contents
			AlignedArray<T> backupCopy(*this);
			resize_zero_shrink(oldNumElements + 1);
			size_t oldNumElementsTripled = 3 * _numEntriesPerArray;
			std::memcpy(this->_p, &(backupCopy[0]), oldNumElementsTripled * sizeof(T));
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
