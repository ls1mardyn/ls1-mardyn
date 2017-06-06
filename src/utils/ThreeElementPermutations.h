/*
 * ThreeElementPermutations.h
 *
 *  Created on: 2 Jun 2017
 *      Author: tchipevn
 */

#ifndef SRC_UTILS_THREEELEMENTPERMUTATIONS_H_
#define SRC_UTILS_THREEELEMENTPERMUTATIONS_H_

#include "mardyn_assert.h"
#include <array>

namespace Permute3Elements {

enum Permutation {
	XYZ,
	XZY,
	YZX,
	YXZ,
	ZXY,
	ZYX
};

template <typename T>
Permutation getPermutationForIncreasingSorting(std::array<T, 3> dims);

template <typename T>
std::array<T, 3> permuteForward(Permutation p, const std::array<T, 3>& in) {
	std::array<T, 3> out;
	switch (p) {
	case XYZ:
		out[0] = in[0]; out[1] = in[1]; out[2] = in[2];
		break;
	case XZY:
		out[0] = in[0]; out[1] = in[2]; out[2] = in[1];
		break;
	case YXZ:
		out[0] = in[1]; out[1] = in[0]; out[2] = in[2];
		break;
	case YZX:
		out[0] = in[1]; out[1] = in[2]; out[2] = in[0];
		break;
	case ZXY:
		out[0] = in[2]; out[1] = in[0]; out[2] = in[1];
		break;
	case ZYX:
		out[0] = in[2]; out[1] = in[1]; out[2] = in[0];
		break;
	}
	return out;
}

template <typename T>
static std::array<T, 3> permuteBackward(Permutation p, const std::array<T, 3>& in) {
	std::array<T, 3> out;

	switch(p) {
	case XYZ:
		out[0] = in[0]; out[1] = in[1]; out[2] = in[2];
		break;
	case XZY:
		out[0] = in[0]; out[2] = in[1]; out[1] = in[2];
		break;
	case YXZ:
		out[0] = in[1]; out[1] = in[0]; out[2] = in[2];
		break;
	case YZX:
		out[0] = in[2]; out[1] = in[0]; out[2] = in[1];
		break;
	case ZXY:
		out[0] = in[1]; out[1] = in[2]; out[2] = in[0];
		break;
	case ZYX:
		out[0] = in[2]; out[1] = in[1]; out[2] = in[0];
		break;
	}
	return out;
}

template <typename T>
Permutation getPermutationForIncreasingSorting(std::array<T, 3> dims) {
	T z_dim = dims[2];
	T y_dim = dims[1];
	T x_dim = dims[0];
	Permutation permutation;

	if (y_dim < x_dim) {
		if (z_dim < x_dim) {
			if (z_dim < y_dim) {
				permutation = ZYX;
			} else {
				permutation = YZX;
			}
		} else {
			permutation = YXZ;
		}
	} else {
		if (z_dim < y_dim) {
			if (z_dim < x_dim) {
				permutation = ZXY;
			} else {
				permutation = XZY;
			}
		} else {
			permutation = XYZ;
		}
	}

#ifndef NDEBUG
	std::array<T, 3> dimsPerm = permuteForward(permutation, dims);
	mardyn_assert(dimsPerm[2] >= dimsPerm[1]);
	mardyn_assert(dimsPerm[1] >= dimsPerm[0]);
#endif

	return permutation;
}

} /* namespace threeElementPermutations */



#endif /* SRC_UTILS_THREEELEMENTPERMUTATIONS_H_ */
