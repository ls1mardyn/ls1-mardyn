/**
 * AlignedAllocator.h
 * Micha Müller
 *
 * with help of:
 * http://coliru.stacked-crooked.com/a/9214894a7c7a921d
 * http://stackoverflow.com/revisions/12942652/2
 */

#ifndef ALIGNEDALLOCATOR_H
#define ALIGNEDALLOCATOR_H

#ifdef __SSE3__
#include <xmmintrin.h>
#endif

#include <limits>
#include <stdlib.h>
//#include <malloc.h>

#define CACHE_LINE_SIZE 64

/**
 * \brief A custom allocator to get aligned memory
 * \details This allocator is intended to be used by std::vector
 * \tparam T The type of the elements this class should allocate memory for
 * \tparam alignment The alignment restriction. Must be a power of 2, should not be 8.
 * \author Micha Mueller
 */
template<typename T, size_t Alignment = CACHE_LINE_SIZE>
struct AlignedAllocator {

	typedef T value_type;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef size_t size_type;

	//dont know what rebind exactly does
	//but makes the allocator compile with second template parameter
	template<class U>
	struct rebind {
		typedef AlignedAllocator<U, Alignment> other;
	};

	/**
	 * \brief Default constructor
	 */
	AlignedAllocator() = default;

	/**
	 * \brief Copy constructor
	 */
	template<class U>
	AlignedAllocator(const AlignedAllocator<U, Alignment>&) {
	}

	/**
	 * \brief Returns maximum possible value of n, with which we can call allocate(n)
	 */
	size_t max_size() const noexcept {
		return (std::numeric_limits < size_t > ::max() - size_t(Alignment)) / sizeof(T);
	}

	/**
	 * \brief Allocate aligned memory for n objects of type T
	 * \return Pointer to the allocated memory
	 */
	T* allocate(std::size_t n) {
		if (n <= max_size()) {
#if defined(_SX)
			T* ptr = static_cast<T*>(malloc(sizeof(T) * n));
#elif defined(__SSE3__) && !defined(__PGI)
			T* ptr = static_cast<T*>(_mm_malloc(sizeof(T) * n, Alignment));
#else
			T* ptr = static_cast<T*>(memalign(Alignment, sizeof(T) * n));
			//T* ptr = static_cast<T*>(aligned_alloc(Alignment, sizeof(T) * n));
			//T* ptr; posix_memalign(&ptr,Alignment, sizeof(T) * n);
#endif
			if (ptr == nullptr) {
				throw std::bad_alloc();
			}
			return ptr;
		}
		throw std::bad_alloc();
	}

	/**
	 * \brief Deallocate memory pointed to by ptr
	 */
	void deallocate(T* ptr, std::size_t /*n*/) {
#if defined(__SSE3__) && !defined(__PGI)
		_mm_free(ptr);
#else
		free(ptr);
#endif
	}

	/**
	 * \brief Construct object of type U at already allocated memory, pointed to by p
	 */
	template<class U, class ...Args>
	void construct(U* p, Args&&... args) noexcept(noexcept(U(std::forward<Args>(args)...))) {
		::new ((void*) p) U(std::forward<Args>(args)...);
	}

	/**
	 * \brief Destroy object pointed to by p, but does not deallocate the memory
	 */
	template<class U>
	void destroy(U* p) noexcept(noexcept(p->~U())) {
		p->~U();
	}
};

template<typename T, size_t TAlignment, typename U, size_t UAlignment>
inline bool operator ==(const AlignedAllocator<T, TAlignment>&,
		const AlignedAllocator<U, UAlignment>&) {
	return TAlignment == UAlignment;
}

template<typename T, size_t TAlignment, typename U, size_t UAlignment>
inline bool operator !=(const AlignedAllocator<T, TAlignment>& a,
		const AlignedAllocator<U, UAlignment>& b) {
	return !(a == b);
}

#endif /*ALIGNEDALLOCATOR_H*/
