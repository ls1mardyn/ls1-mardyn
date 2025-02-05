/*
 * RealAccumVecTest.cpp
 *
 *  Created on: 10 Jan 2018
 *      Author: tchipevn
 */

#include "RealAccumVecTest.h"
#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"
#include "particleContainer/adapter/vectorization/SIMD_DEFINITIONS.h"
#include "particleContainer/adapter/vectorization/MaskVec.h"
#include "utils/AlignedArray.h"

TEST_SUITE_REGISTRATION(RealAccumVecTest);

RealAccumVecTest::RealAccumVecTest() {
}

RealAccumVecTest::~RealAccumVecTest() {}

void RealAccumVecTest::setUp() {}
void RealAccumVecTest::tearDown() {}

void RealAccumVecTest::testConvert() {
	AlignedArray<vcp_real_calc> calcArray(VCP_VEC_SIZE);
	for (unsigned i = 0; i < VCP_VEC_SIZE; ++i) {
		calcArray.appendValue(static_cast<vcp_real_calc>(i), i);
	}
	RealCalcVec calc = RealCalcVec::aligned_load(calcArray);

	RealAccumVec accum;
	accum = RealAccumVec::zero();
	accum = RealAccumVec::convertCalcToAccum(calc);

	AlignedArray<vcp_real_accum> accumArray(VCP_VEC_SIZE);
	vcp_real_accum * const accumArrayData = accumArray;

	accum.aligned_store(accumArrayData);

	for (unsigned i = 0; i < VCP_VEC_SIZE; ++i) {
		ASSERT_DOUBLES_EQUAL(static_cast<vcp_real_accum>(i),accumArray[i],1e-16);
	}

	RealAccumVec accum2 = RealAccumVec::aligned_load(accumArrayData);

	AlignedArray<vcp_real_accum> accumArray2(VCP_VEC_SIZE);
	vcp_real_accum * const accumArrayData2 = accumArray2;
	accum2.aligned_store(accumArrayData2);

	for (unsigned i = 0; i < VCP_VEC_SIZE; ++i) {
		ASSERT_DOUBLES_EQUAL(accumArray[i],accumArray2[i],1e-16);
	}


	test_log -> info() << " ok" << std::endl;
}

