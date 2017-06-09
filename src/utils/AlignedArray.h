/**
 * \file
 * \brief AlignedArray.h
 */

#ifndef ALIGNEDARRAY_H_
#define ALIGNEDARRAY_H_

#ifdef __SSE3__
#include <xmmintrin.h>
#endif
#include <malloc.h>
#include <new>
#include <cstring>
#include <vector>
#include "utils/mardyn_assert.h"
#include "AlignedAllocator.h"

#define CACHE_LINE_SIZE 64

// TODO: managing this class is becoming too complicated
// mainly because of the need to remember what happens when,
// since our functions don't follow standard naming rules
// and conventions stemming from std::vector.
//
// Switch to std::vector with a custom allocator
// so that it is at least clear what happens when.

/**
 * \brief An aligned array.
 * \details Has pointer to T semantics. Internal size is rounded up to fill up full cache-lines.
 * \tparam T The type of the array elements.
 * \tparam alignment The alignment restriction. Must be a power of 2, should not be 8.
 * \author Johannes Heckl, Nikola Tchipev, Micha Mueller
 */
template<class T, size_t alignment = CACHE_LINE_SIZE>
class AlignedArray {
public:
	/**
	 * \brief Construct an empty array.
	 */
	AlignedArray() :
			_vec(0) {
	}

	/**
	 * \brief Construct an array of n elements.
	 */
	AlignedArray(size_t n) :
			_vec(n) {
	}

	/**
	 * \brief Construct a copy of another AlignedArray.
	 */
	AlignedArray(const AlignedArray & a) :
			_vec(a._vec) {
	}

	/**
	 * \brief Assign a copy of another AlignedArray.
	 */
	AlignedArray & operator=(const AlignedArray & a) {
		_vec = a._vec;
		return *this;
	}

	/**
	 * \brief Free the array.
	 */
	virtual ~AlignedArray() {
	}

	void appendValue(T v, size_t oldNumElements) {
		mardyn_assert(oldNumElements <= _vec.size());
		if (oldNumElements >= _vec.size()) {
			// shit, we need to resize, but also keep contents
			_vec.resize(_round_up(oldNumElements + 1));
		}
		_vec[oldNumElements] = v;
	}

	virtual size_t resize_zero_shrink(size_t exact_size, bool zero_rest_of_CL =
			false, bool allow_shrink = false) {
		size_t size_rounded_up = _round_up(exact_size);

		bool need_resize = size_rounded_up > _vec.size()
				or (allow_shrink and size_rounded_up < _vec.size());

		if (need_resize) {
			_vec.resize(size_rounded_up);
		}
		// we might still need to zero the rest of the Cache Line
		if (zero_rest_of_CL and size_rounded_up > 0) {
			std::memset(_vec.data() + exact_size, 0,
					size_rounded_up - exact_size);
		}

		mardyn_assert(size_rounded_up <= _vec.size());
		return _vec.size();
	}

	/**
	 * \brief Reallocate the array. All content may be lost.
	 */
	virtual void resize(size_t n) {
		_vec.resize(_round_up(n));
	}

	virtual void zero(size_t start_idx) {
		if (_vec.size() > 0) {
			size_t num_to_zero = this->_round_up(start_idx) - start_idx;
			std::memset(_vec.data(), 0, num_to_zero * sizeof(T));
		}
	}

	/**
	 * \brief Return current size in terms of elements
	 */
	inline size_t get_size() const {
		return _vec.size();
	}

	/**
	 * \brief Implicit conversion into pointer to T.
	 */
	operator T*() {
		return _vec.data();
	}

	operator const T*() const {
		return _vec.data();
	}

	/**
	 * \brief Return amount of allocated storage + .
	 */
	size_t get_dynamic_memory() const {
		return _vec.size() * sizeof(T);
	}

	static size_t _round_up(size_t n) {
		size_t ret = 0;
		switch (sizeof(T)) {
		case 1:
			ret = (n + 63) & ~0x3F;
			break;
		case 2:
			ret = (n + 31) & ~0x1F;
			break;
		case 4:
			ret = (n + 15) & ~0x0F;
			break;
		case 8:
			ret = (n + 7) & ~0x07;
			break;
		default:
			mardyn_assert(false);
		}
		return ret;
	}

protected:

	std::vector<T, AlignedAllocator<T, alignment>> _vec;

};

#endif

