/*
 *  Created on: 30 May 2017
 *      Author: Micha Mueller
 */

#include "RealCalcVecTest.h"
#include "particleContainer/adapter/vectorization/SIMD_DEFINITIONS.h"
#include "utils/AlignedArray.h"

#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER
	TEST_SUITE_REGISTRATION(RealCalcVecTest);
#else
	#pragma message "Compilation info: The unit tests for the fast reciproce methods are only executed with AVX2, KNL and AVX512"
#endif


RealCalcVecTest::RealCalcVecTest() {
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or VCP_VEC_TYPE == VCP_VEC_KNL or VCP_VEC_TYPE == VCP_VEC_KNL_GATHER
	test_log->info() << "Testing RealCalcVec fast reciprocal with "
	#if VCP_VEC_TYPE==VCP_VEC_AVX2
		<< "AVX2 intrinsics." << std::endl;
	#elif VCP_VEC_TYPE==VCP_VEC_KNL
		<< "KNL_MASK intrinsics." << std::endl;
	#elif VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
		<< "KNL_G_S intrinsics." << std::endl;
	#elif VCP_VEC_TYPE==VCP_VEC_AVX512F
		<< "SKX_MASK intrinsics." << std::endl;
	#elif VCP_VEC_TYPE==VCP_VEC_AVX512F_GATHER
		<< "SKX_G_S intrinsics." << std::endl;
	#endif
#endif
}

RealCalcVecTest::~RealCalcVecTest() {}

void RealCalcVecTest::setUp() {}

void RealCalcVecTest::tearDown() {}

void RealCalcVecTest::testFastReciprocalMask() {
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER

	AlignedArray<vcp_real_calc> storageLegacy(VCP_VEC_SIZE);
	storageLegacy.appendValue(1.0, 0);
	storageLegacy.appendValue(-1e-6, 1);
	storageLegacy.appendValue(1234999.56, 2);
	storageLegacy.appendValue(0.0, 3);

	AlignedArray<vcp_real_calc> storageFast(VCP_VEC_SIZE);
	storageFast.appendValue(1.0, 0);
	storageFast.appendValue(-1e-6, 1);
	storageFast.appendValue(1234999.56, 2);
	storageFast.appendValue(0.0, 3);

	#if VCP_VEC_TYPE == VCP_VEC_AVX2
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			vcp::MaskCalcVec mask = _mm256_set_epi32(0, 0, 0, 0, 0, ~0, ~0, ~0);
		#else /* VCP_DPDP */
			vcp::MaskCalcVec mask = _mm256_set_epi64x(0, ~0, ~0, ~0);
		#endif
	#else
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			vcp::MaskCalcVec mask = 0x0007;
		#else /* VCP_DPDP */
			vcp::MaskCalcVec mask = 0x07;
		#endif
	#endif

	RealCalcVec test = RealCalcVec::aligned_load(storageLegacy);
	RealCalcVec result_unmasked = RealCalcVec::set1(1.0) / test;
	RealCalcVec result = RealCalcVec::apply_mask(result_unmasked, mask);
	result.aligned_store(storageLegacy);

	test = RealCalcVec::aligned_load(storageFast);
	result = RealCalcVec::fastReciprocal_mask(test, mask);
	result.aligned_store(storageFast);

	for(size_t i = 0; i < VCP_VEC_SIZE; ++i) {
		std::stringstream str;
		str << "Value at index i = " << i << std::endl;
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), storageLegacy[i], storageFast[i], 1e-6);
		#else /* VCP_DPDP */
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), storageLegacy[i], storageFast[i], 1e-8);
		#endif
	}
#endif
}

void RealCalcVecTest::testFastReciprocSqrtMask() {
#if VCP_VEC_TYPE == VCP_VEC_AVX2 or \
	VCP_VEC_TYPE == VCP_VEC_KNL or \
	VCP_VEC_TYPE == VCP_VEC_KNL_GATHER or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F or \
	VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER

	AlignedArray<vcp_real_calc> storageLegacy(VCP_VEC_SIZE);
	storageLegacy.appendValue(1.0, 0);
	storageLegacy.appendValue(1e-6, 1);
	storageLegacy.appendValue(1234999.56, 2);
	storageLegacy.appendValue(0.0, 3);

	AlignedArray<vcp_real_calc> storageFast(VCP_VEC_SIZE);
	storageFast.appendValue(1.0, 0);
	storageFast.appendValue(1e-6, 1);
	storageFast.appendValue(1234999.56, 2);
	storageFast.appendValue(0.0, 3);

	#if VCP_VEC_TYPE == VCP_VEC_AVX2
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			vcp::MaskCalcVec mask = _mm256_set_epi32(0, 0, 0, 0, 0, ~0, ~0, ~0);
		#else /* VCP_DPDP */
			vcp::MaskCalcVec mask = _mm256_set_epi64x(0, ~0, ~0, ~0);
		#endif
	#else
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			vcp::MaskCalcVec mask = 0x0007;
		#else /* VCP_DPDP */
			vcp::MaskCalcVec mask = 0x07;
		#endif
	#endif

	RealCalcVec test = RealCalcVec::aligned_load(storageLegacy);
	RealCalcVec result_unmasked = RealCalcVec::set1(1.0) / RealCalcVec::sqrt(test);
	RealCalcVec result = RealCalcVec::apply_mask(result_unmasked, mask);
	result.aligned_store(storageLegacy);

	test = RealCalcVec::aligned_load(storageFast);
	result = RealCalcVec::fastReciprocSqrt_mask(test, mask);
	result.aligned_store(storageFast);

	for(size_t i = 0; i < VCP_VEC_SIZE; ++i) {
		std::stringstream str;
		str << "Value at index i = " << i << std::endl;
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), storageLegacy[i], storageFast[i], 1e-4);
		#else /* VCP_DPDP */
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), storageLegacy[i], storageFast[i], 1e-8);
		#endif
	}
#endif
}
