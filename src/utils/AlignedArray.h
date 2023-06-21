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
//
// Regarding prefetching on Xeon Phi, consider the following text, coming from here:
// https://tianyuliukingcrimson.wordpress.com/2015/07/01/prefetch-on-intel-mic-coprocessor-and-xeon-cpu/
	//Prefetch instruction
	//
	//Let’s take a look at two orthogonal concepts first:
	//
	//    non-temporal hint (NTA) — informs that data will be used only once in the future and causes them to be evicted from the cache after the first use (most recently used data to be evicted).
	//    exclusive hint (E) — renders the cache line on the current core in the “exclusive” state, where the cache lines on other cores are invalidated.
	//
	//The combination of temporality, exclusiveness, and locality (L1 or L2) together yields 8 types of instructions supported by the present-day Knights Corner MIC. They specify how the data are expected to be uniquely handled in the cache, enumerated below.
	//instruction 	hint 	purpose
	//vprefetchnta 	_MM_HINT_NTA 	loads data to L1 and L2 cache, marks it as NTA
	//vprefetch0 	_MM_HINT_T0 	loads data to L1 and L2 cache
	//vprefetch1 	_MM_HINT_T1 	loads data to L2 cache only
	//vprefetch2 	_MM_HINT_T2 	loads data to L2 cache only, marks it as NTA This mnemonic is counter-intuitive as there is not NTA in it
	//vprefetchenta 	_MM_HINT_ENTA 	exclusive version of vprefetchnta
	//vprefetche0 	_MM_HINT_ET0 	exclusive version of vprefetch0
	//vprefetche1 	_MM_HINT_ET1 	exclusive version of vprefetch1
	//vprefetche2 	_MM_HINT_ET2 	exclusive version of vprefetch2
	//
	//Note L2 cache of the MIC is inclusive in the sense that it has a copy of all the data in L1.
	//
	//There are two ways of implementing prefetch in C — intrinsic and assembly.
	//
	//// method 1: intrinsic
	//_mm_prefetch((const char*)addr, hint);
	//
	//// method 2: assembly
	//asm volatile ("prefetch_inst [%0]"::"m"(addr));
	//
	//Here addr is the address of the byte starting from which to prefetch, prefetch_inst is the prefetch instructions listed above, and hint is the parameter for the compiler intrinsic. We would like to emphasize again that _MM_HINT_T2 and _MM_HINT_ET2 are counter-intuitive. In fact they are misnomers as both are non-temporary. They should have been named as _MM_HINT_NTA2 and _MM_HINT_ENTA2 by Intel.

// Prefetching on Xeons features much fewer hints, apparently! (same link):
	//prefetchnta 	_MM_HINT_NTA 	loads data to L2 and L3 cache, marks as NTA
	//prefetcht0 	_MM_HINT_T0 	loads data to L2 and L3 cache
	//prefetcht1 	_MM_HINT_T1 	equivalent to prefetch0
	//prefetcht2 	_MM_HINT_T2 	equivalent to prefetch0

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

#if defined(__SSE3__) or defined(__MIC__)
	virtual void prefetch(int /*hint*/ = 1, int n = -1) const {
		mardyn_assert(n >= -2);

		size_t endPrefetch;
		const size_t stride = _round_up(1);

		switch(n) {
		case -1:
			// prefetch all up to capacity()
			endPrefetch = _vec.capacity();
			break;
		case -2:
			// prefetch all up to size()
			endPrefetch = _vec.size();
			break;
		default:
			// prefetch only first n elements
			endPrefetch = static_cast<size_t>(n);
		}

		for (size_t i = 0; i < endPrefetch; i+= stride) {
			const T & val = _vec[i];
			const T * valP = &val;
#if defined(__MIC__)
			_mm_prefetch((const char*)valP, 2);
#else
			_mm_prefetch((const char*)valP, _MM_HINT_T1);
#endif
		}
	}
#else
	virtual void prefetch(int /*hint = 1*/, int /*n = -1*/) const {}
#endif

	virtual void increaseStorage(size_t oldNumElements, size_t additionalElements) {
		mardyn_assert(oldNumElements <= _vec.capacity());

		size_t newNumElements = oldNumElements + additionalElements;

		if (newNumElements <= _vec.capacity()) {
			// no need to resize
			return;
		}

		// we need to resize, but also keep contents
		_vec.reserve(_round_up(newNumElements));
		_vec.resize(_vec.capacity());
	}

	void appendValue(T v, size_t oldNumElements) {
		increaseStorage(oldNumElements, 1);

		_vec[oldNumElements] = v;
	}

	virtual size_t resize_zero_shrink(size_t exact_size, bool zero_rest_of_CL =
			false, bool allow_shrink = false) {
		size_t size_rounded_up = _round_up(exact_size);

		bool need_resize = size_rounded_up > _vec.size()
				or (allow_shrink and size_rounded_up < _vec.size());

		if (need_resize) {
			_vec.reserve(size_rounded_up);
			_vec.resize(_vec.capacity());
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
		_vec.reserve(_round_up(n));
		_vec.resize(_vec.capacity());
	}

	virtual void zero(size_t start_idx = 0) {
		if (_vec.size() > 0 and start_idx < _vec.capacity()) {
			size_t num_to_zero = _vec.capacity() - start_idx;
			std::memset(_vec.data() + start_idx, 0, num_to_zero * sizeof(T));
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
		return _vec.capacity() * sizeof(T);
	}

	static size_t _round_up(size_t n) {
		unsigned long j = alignment / sizeof(T) - 1;
		unsigned long ret = (n + j) & ~j;
		return ret;
	}

protected:

	std::vector<T, AlignedAllocator<T, alignment>> _vec;

};

#endif

