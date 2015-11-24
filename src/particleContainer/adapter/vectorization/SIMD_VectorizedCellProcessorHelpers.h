/***********************************************************************************//**
 *
 * \file SIMD_TYPES.h
 *
 * \brief Defines the length of the vectors and the corresponding functions
 *
 * \author Steffen Seckler
 *
 **************************************************************************************/
#pragma once
#ifndef  SIMD_VECTORIZEDCELLPROCESSORHELPERS_H
#define  SIMD_VECTORIZEDCELLPROCESSORHELPERS_H

#include "SIMD_TYPES.h"
#include "utils/AlignedArray.h"
typedef AlignedArray<double> DoubleArray;

static inline __attribute__((always_inline))
void unpackEps24Sig2(vcp_double_vec& eps_24, vcp_double_vec& sig2, const DoubleArray& eps_sigI,
		const size_t* const id_j){
#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	eps_24 = eps_sigI[2 * id_j[0]];
	sig2 = eps_sigI[2 * id_j[0] + 1];

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	const vcp_double_vec e0s0 = vcp_simd_load(eps_sigI + 2 * id_j[0]);
	const vcp_double_vec e1s1 = vcp_simd_load(eps_sigI + 2 * id_j[1]);
	eps_24 = vcp_simd_unpacklo(e0s0, e1s1);
	sig2 = vcp_simd_unpackhi(e0s0, e1s1);

#elif VCP_VEC_TYPE==VCP_VEC_AVX //avx
	static const __m256i memoryMask_first_second = _mm256_set_epi32(0, 0, 0, 0, 1<<31, 0, 1<<31, 0);
	const vcp_double_vec e0s0 = vcp_simd_maskload(eps_sigI + 2 * id_j[0], memoryMask_first_second);
	const vcp_double_vec e1s1 = vcp_simd_maskload(eps_sigI + 2 * id_j[1], memoryMask_first_second);
	const vcp_double_vec e2s2 = vcp_simd_maskload(eps_sigI + 2 * id_j[2], memoryMask_first_second);
	const vcp_double_vec e3s3 = vcp_simd_maskload(eps_sigI + 2 * id_j[3], memoryMask_first_second);

	const vcp_double_vec e0e1 = vcp_simd_unpacklo(e0s0, e1s1);
	const vcp_double_vec s0s1 = vcp_simd_unpackhi(e0s0, e1s1);
	const vcp_double_vec e2e3 = vcp_simd_unpacklo(e2s2, e3s3);
	const vcp_double_vec s2s3 = vcp_simd_unpackhi(e2s2, e3s3);

	eps_24 = _mm256_permute2f128_pd(e0e1, e2e3, 1<<5);
	sig2 = _mm256_permute2f128_pd(s0s1, s2s3, 1<<5);
#endif
}

static inline __attribute__((always_inline))
void unpackShift6(vcp_double_vec& shift6, const DoubleArray& shift6I,
		const size_t* const id_j){
#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	shift6 = shift6I[id_j[0]];

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	const vcp_double_vec sh1 = _mm_load_sd(shift6I + id_j[0]);
	const vcp_double_vec sh2 = _mm_load_sd(shift6I + id_j[1]);
	shift6 = vcp_simd_unpacklo(sh1, sh2);

#elif VCP_VEC_TYPE==VCP_VEC_AVX //avx
	static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
	const vcp_double_vec sh0 = vcp_simd_maskload(shift6I + id_j[0], memoryMask_first);
	const vcp_double_vec sh1 = vcp_simd_maskload(shift6I + id_j[1], memoryMask_first);
	const vcp_double_vec sh2 = vcp_simd_maskload(shift6I + id_j[2], memoryMask_first);
	const vcp_double_vec sh3 = vcp_simd_maskload(shift6I + id_j[3], memoryMask_first);

	const vcp_double_vec sh0sh1 = vcp_simd_unpacklo(sh0, sh1);
	const vcp_double_vec sh2sh3 = vcp_simd_unpacklo(sh2, sh3);

	shift6 = _mm256_permute2f128_pd(sh0sh1, sh2sh3, 1<<5);

#endif
}




























#endif /* SIMD_VECTORIZEDCELLPROCESSORHELPERS_H */
