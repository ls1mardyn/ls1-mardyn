#ifndef UNORDERED_VECTOR_H
#define UNORDERED_VECTOR_H

#include <vector>
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
//! @param pos iterator pointing to the position at which we want to erase (in)
//! @return iterator to "next" element
template<typename T, typename A>
void fastRemove(std::vector<T, A>& v, typename std::vector<T, A>::iterator& pos) {
	using std::swap;

	// assumption: if the vector is empty, then the method is not called at all
	// i.e. v is not empty()
	mardyn_assert(not v.empty());

	// assert that the iterator points inside this vector
	mardyn_assert(pos >= v.begin());
	mardyn_assert(pos < v.end());

	// either std::swap, or a user-written function in the namespace
	swap(*pos, v.back());
	// note: swapping is necessary (as opposed to only copying )

	v.pop_back(); // end() will change and may become equal to pos, which is intended
}

} /* namespace UnorderedVector */

#endif /* UNORDERED_VECTOR_H */
