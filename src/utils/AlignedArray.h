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
#include <cassert>

#define CACHE_LINE_SIZE 64

/**
 * \brief An aligned array.
 * \details Has pointer to T semantics.
 * \tparam T The type of the array elements.
 * \tparam alignment The alignment restriction. Must be a power of 2, should not be 8.
 * \author Johannes Heckl, Nikola Tchipev
 */
template<class T, size_t alignment = CACHE_LINE_SIZE>
class AlignedArray {
public:
	/**
	 * \brief Construct an empty array.
	 */
	AlignedArray() :
			_n(0), _p(nullptr) {
	}

	/**
	 * \brief Construct an array of n elements.
	 */
	AlignedArray(size_t n) :
			_n(n), _p(_allocate(n)) {
		if (!_p)
			throw std::bad_alloc();
	}

	/**
	 * \brief Construct a copy of another AlignedArray.
	 */
	AlignedArray(const AlignedArray & a) :
			_n(a._n), _p(_allocate(a._n)) {
		if (!_p)
			throw std::bad_alloc();
		_assign(a);
	}

	/**
	 * \brief Assign a copy of another AlignedArray.
	 */
	AlignedArray & operator=(const AlignedArray & a) {
		resize(a._n);
		_assign(a);
		return *this;
	}

	/**
	 * \brief Free the array.
	 */
	virtual ~AlignedArray() {
		_free();
	}

	virtual size_t resize_zero_shrink(size_t exact_size, bool zero_rest_of_CL = false, bool allow_shrink = false) {
		size_t size_rounded_up = _round_up(exact_size);

		bool need_resize = size_rounded_up > _n or (allow_shrink and size_rounded_up < _n);

		if (need_resize) {
			resize(size_rounded_up);
			// resize zero-s all
		} else {
			// we didn't resize, but we might still need to zero the rest of the Cache Line
			if (zero_rest_of_CL and size_rounded_up > 0) {
				std::memset(_p + exact_size, 0, size_rounded_up - exact_size);
			}
		}

		assert(size_rounded_up == _n);
		return _n;
	}

	/**
	 * \brief Reallocate the array. All content may be lost.
	 */
	virtual void resize(size_t n) {
		if (n == _n)
			return;
		_free();
		_p = nullptr;
		_p = _allocate(n);
		if (_p == nullptr)
			throw std::bad_alloc();
		_n = n;
	}

	inline size_t get_size() const {
		return this->_n;
	}

	/**
	 * \brief Implicit conversion into pointer to T.
	 */
	operator T*() const {
		return _p;
	}

	/**
	 * \brief Return amount of allocated storage + .
	 */
	size_t get_dynamic_memory() const {
		return _n * sizeof(T);
	}

	static size_t _round_up(size_t n) {
		// NOTE: this function is never called
		// there are specializations provided after the class, as
		// this function may be called in performance intensive places
		assert(false);

		size_t multiple = CACHE_LINE_SIZE / sizeof(T);
		assert(multiple * sizeof(T) == CACHE_LINE_SIZE);
		return ((n + multiple - 1) / multiple) * multiple;
	}

protected:
	void _assign(T * p) const {
		std::memcpy(_p, p, _n * sizeof(T));
	}

	static T* _allocate(size_t elements) {

#if defined(__SSE3__) && ! defined(__PGI)
		T* ptr = static_cast<T*>(_mm_malloc(sizeof(T) * elements, alignment));
#else
		T* ptr = static_cast<T*>(memalign(alignment, sizeof(T) * elements));
#endif

		std::memset(ptr, 0, elements * sizeof(T));
		return ptr;

	}

	void _free()
	{
#if defined(__SSE3__) && ! defined(__PGI)
		_mm_free(_p);
#else
		free(_p);
#endif
	}

	size_t _n;
	T * _p;
};

template<>
inline size_t AlignedArray<double>	::_round_up(size_t n) {
	return (n + 7) & ~0x07;
}

template<>
inline size_t AlignedArray<size_t>	::_round_up(size_t n) {
	return (n + 7) & ~0x07;
}

template<>
inline size_t AlignedArray<float>	::_round_up(size_t n) {
	return (n + 15) & ~0x0F;
}

template<>
inline size_t AlignedArray<int>		::_round_up(size_t n) {
	return (n + 15) & ~0x0F;
}



#endif
