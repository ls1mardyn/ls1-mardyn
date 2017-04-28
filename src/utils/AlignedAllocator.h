/**
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

#define CACHE_LINE_SIZE 64

/**
 * Custom allocator for allocating aligned memory.
 */
template <typename T, size_t Alignment = CACHE_LINE_SIZE>
struct AlignedAllocator {

        using value_type = T;
        using pointer = T*;

        //dont know what rebind is for
        //but makes the allocator compile with second template parameter
        template <class U>
        struct rebind { typedef AlignedAllocator<U, Alignment> other; };

        //default constructor
        AlignedAllocator() = default;

        //copy constructor
        template <class U>
        AlignedAllocator(const AlignedAllocator<U, Alignment>&) {
        }

        T* allocate(std::size_t n) {
                if (n <= std::numeric_limits<std::size_t>::max() / sizeof(T)) {
                        #if defined(__SSE3__) && !defined(__PGI)
                        T* ptr = static_cast<T*>(_mm_malloc(sizeof(T) * n, Alignment));
                        #else
                        T* ptr = static_cast<T*>(memalign(Alignment, sizeof(T) * n));
                        #endif
                        if(ptr == nullptr) {
                          throw std::bad_alloc();
                        }
                        return ptr;
                }
                throw std::bad_alloc();
        }

        void deallocate(T* ptr, std::size_t n) {
                #if defined(__SSE3__) && !defined(__PGI)
                _mm_free(ptr);
                #else
                free(ptr);
                #endif
        }
};

template <typename T, size_t TAlignment, typename U, size_t UAlignment>
inline bool operator == (const AlignedAllocator<T, TAlignment>&, const AlignedAllocator<U, UAlignment>&) {
        return TAlignment == UAlignment;
}

template <typename T, size_t TAlignment, typename U, size_t UAlignment>
inline bool operator != (const AlignedAllocator<T, TAlignment>& a, const AlignedAllocator<U, UAlignment>& b) {
        return !(a == b);
}

#endif /*ALIGNEDALLOCATOR_H*/
