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

	T* xBegin() { return _numEntriesPerArray > 0 ? this->_vec.data() + (0 * _numEntriesPerArray) : nullptr; }
	T* yBegin() { return _numEntriesPerArray > 0 ? this->_vec.data() + (1 * _numEntriesPerArray) : nullptr; }
	T* zBegin() { return _numEntriesPerArray > 0 ? this->_vec.data() + (2 * _numEntriesPerArray) : nullptr; }
	T* xBegin() const { return _numEntriesPerArray > 0 ? this->_vec_ptr->data() + (0 * _numEntriesPerArray) : nullptr; }
	T* yBegin() const { return _numEntriesPerArray > 0 ? this->_vec_ptr->data() + (1 * _numEntriesPerArray) : nullptr; }
	T* zBegin() const { return _numEntriesPerArray > 0 ? this->_vec_ptr->data() + (2 * _numEntriesPerArray) : nullptr; }

	T& x(size_t i) { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 0 * _numEntriesPerArray]; }
	T& y(size_t i) { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 1 * _numEntriesPerArray]; }
	T& z(size_t i) { mardyn_assert(i < _numEntriesPerArray); return this->_vec[i + 2 * _numEntriesPerArray]; }
	T& x(size_t i) const { mardyn_assert(i < _numEntriesPerArray); return (*this->_vec_ptr)[i + 0 * _numEntriesPerArray]; }
	T& y(size_t i) const { mardyn_assert(i < _numEntriesPerArray); return (*this->_vec_ptr)[i + 1 * _numEntriesPerArray]; }
	T& z(size_t i) const { mardyn_assert(i < _numEntriesPerArray); return (*this->_vec_ptr)[i + 2 * _numEntriesPerArray]; }

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

	void appendValueTriplet(T v0, T v1, T v2, size_t oldNumElements) {
		mardyn_assert(oldNumElements <= _numEntriesPerArray);
		if (oldNumElements < _numEntriesPerArray) {
			// no need to resize, baby
		} else {
			// shit, we need to resize

			// do we need to keep contents?
			if(oldNumElements == 0) {
				// No
				resize_zero_shrink(oldNumElements + 1);
			} else {
				// Yes
				AlignedArray<T> backupCopy(*this);
				resize_zero_shrink(oldNumElements + 1);
				std::memcpy(xBegin(), &(backupCopy[0*oldNumElements]), oldNumElements * sizeof(T));
				std::memcpy(yBegin(), &(backupCopy[1*oldNumElements]), oldNumElements * sizeof(T));
				std::memcpy(zBegin(), &(backupCopy[2*oldNumElements]), oldNumElements * sizeof(T));
			}
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

