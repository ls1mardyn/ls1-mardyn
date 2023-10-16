/*
 * PermutationTest.cpp
 *
 *  Created on: 2 Jun 2017
 *      Author: tchipevn
 */

#include "PermutationTest.h"
#include "../ThreeElementPermutations.h"

TEST_SUITE_REGISTRATION(PermutationTest);

void PermutationTest::testPermutations() {
	using std::array;
	using namespace Permute3Elements;

	Log::global_log->info() << "Testing ThreeElementPermutations." << std::endl ;

	std::array<int, 3> a_012 = {0, 1, 2};
	std::array<int, 3> a_021 = {0, 2, 1};
	std::array<int, 3> a_102 = {1, 0, 2};
	std::array<int, 3> a_120 = {1, 2, 0};
	std::array<int, 3> a_201 = {2, 0, 1};
	std::array<int, 3> a_210 = {2, 1, 0};

	Permutation p_012 = getPermutationForIncreasingSorting(a_012);
	Permutation p_021 = getPermutationForIncreasingSorting(a_021);
	Permutation p_102 = getPermutationForIncreasingSorting(a_102);
	Permutation p_120 = getPermutationForIncreasingSorting(a_120);
	Permutation p_201 = getPermutationForIncreasingSorting(a_201);
	Permutation p_210 = getPermutationForIncreasingSorting(a_210);

	ASSERT_EQUAL(p_012, XYZ);
	ASSERT_EQUAL(p_021, XZY);
	ASSERT_EQUAL(p_102, YXZ);
	ASSERT_EQUAL(p_120, ZXY); // note weird
	ASSERT_EQUAL(p_201, YZX); // note weird
	ASSERT_EQUAL(p_210, ZYX);


	// test the forward and backward permutations

	std::array<int, 3> v = {5, 6, -1};
	std::array<int, 3> v1, v2;

	v1 = permuteForward (p_012, v);
	v2 = permuteBackward(p_012, v1);
	for (int d = 0; d < 3; ++d) {
		ASSERT_EQUAL(v2[d], v[d]);
	}

	v1 = permuteForward (p_021, v);
	v2 = permuteBackward(p_021, v1);
	for (int d = 0; d < 3; ++d) {
		ASSERT_EQUAL(v2[d], v[d]);
	}

	v1 = permuteForward (p_102, v);
	v2 = permuteBackward(p_102, v1);
	for (int d = 0; d < 3; ++d) {
		ASSERT_EQUAL(v2[d], v[d]);
	}

	v1 = permuteForward (p_120, v);
	v2 = permuteBackward(p_120, v1);
	for (int d = 0; d < 3; ++d) {
		ASSERT_EQUAL(v2[d], v[d]);
	}

	v1 = permuteForward (p_201, v);
	v2 = permuteBackward(p_201, v1);
	for (int d = 0; d < 3; ++d) {
		ASSERT_EQUAL(v2[d], v[d]);
	}

	v1 = permuteForward (p_210, v);
	v2 = permuteBackward(p_210, v1);
	for (int d = 0; d < 3; ++d) {
		ASSERT_EQUAL(v2[d], v[d]);
	}

}
