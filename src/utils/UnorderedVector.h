#ifndef UNORDERED_VECTOR_H
#define UNORDERED_VECTOR_H

#include <vector>
#include <type_traits>
#include "utils/mardyn_assert.h"

/**
 * \brief Leveraging performance (and especially OpenMP)
 * requires switching from std::lists to std::vectors.
 */
namespace UnorderedVector {

//! @brief In the cases, where the order of insertion is irrelevant,
//!        we can optimize removal from the middle of the list to avoid
//!		   reallocating all subsequent elements by swapping the element
//! 	   to be removed with the last one and calling pop_back()
//! @param v the vector in question (inout)
//! @param index of the position at which we want to erase (in)
template<typename T, typename A, typename UIntType>
void fastRemove(std::vector<T, A>& v, UIntType index) {
	using std::swap;

	mardyn_assert(std::is_integral<UIntType>::value);

	// assumption: if the vector is empty, then the method is not called at all
	// i.e. v is not empty()
	mardyn_assert(not v.empty());

	// assert that the iterator points inside this vector
	mardyn_assert(index >= 0);
	mardyn_assert(index < v.size());

	// either std::swap, or a user-written function in the namespace
	swap(v[index], v.back());
	// note: swapping is necessary (as opposed to only copying )

	v.pop_back(); // end() will change and may become equal to pos, which is intended
}

} /* namespace UnorderedVector */

#endif /* UNORDERED_VECTOR_H */
