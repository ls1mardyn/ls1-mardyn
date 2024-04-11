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

/**
 * unpacks eps_24 and sig2 from the eps_sigI array according to the index array id_j (for mic+avx2: use gather)
 * @param eps_24 vector in which eps_24 is saved
 * @param sig2 vector in which sig2 is saved
 * @param eps_sigI initial eps_sig array
 * @param id_j array of displacements
 * @param offset offset of the id_j array
 */
template <class MaskGatherChooser>
static vcp_inline
void unpackEps24Sig2(RealCalcVec& eps_24, RealCalcVec& sig2, const AlignedArray<vcp_real_calc>& eps_sigI,
					 const vcp_center_id_t* const id_j, const vcp_center_id_t& offset, const vcp_lookupOrMask_vec& lookupORforceMask __attribute__((unused))) {

#if VCP_VEC_TYPE != VCP_VEC_KNL_GATHER and VCP_VEC_TYPE != VCP_VEC_AVX512F_GATHER
	const vcp_center_id_t* id_j_shifted = id_j + offset;//this is the pointer, to where the stuff is stored.
#endif

#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	eps_24 = eps_sigI[2 * id_j_shifted[0]];
	sig2 = eps_sigI[2 * id_j_shifted[0] + 1];

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3

	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		const RealCalcVec e0s0 = _mm_loadl_pi(RealCalcVec::zero(), (const __m64 *)(eps_sigI + 2 * id_j_shifted[0]));
		const RealCalcVec e1s1 = _mm_loadl_pi(RealCalcVec::zero(), (const __m64 *)(eps_sigI + 2 * id_j_shifted[1]));
		const RealCalcVec e2s2 = _mm_loadl_pi(RealCalcVec::zero(), (const __m64 *)(eps_sigI + 2 * id_j_shifted[2]));
		const RealCalcVec e3s3 = _mm_loadl_pi(RealCalcVec::zero(), (const __m64 *)(eps_sigI + 2 * id_j_shifted[3]));

		const RealCalcVec e0e1s0s1 = RealCalcVec::unpack_lo(e0s0, e1s1);
		const RealCalcVec e2e3s2s3 = RealCalcVec::unpack_lo(e2s2, e3s3);

		eps_24 = _mm_castpd_ps(_mm_unpacklo_pd(_mm_castps_pd(e0e1s0s1), _mm_castps_pd(e2e3s2s3)));
		sig2 = _mm_castpd_ps(_mm_unpackhi_pd(_mm_castps_pd(e0e1s0s1), _mm_castps_pd(e2e3s2s3)));

	#else /* VCP_DPDP */
		const RealCalcVec e0s0 = RealCalcVec::aligned_load(eps_sigI + 2 * id_j_shifted[0]);
		const RealCalcVec e1s1 = RealCalcVec::aligned_load(eps_sigI + 2 * id_j_shifted[1]);
		eps_24 = RealCalcVec::unpack_lo(e0s0, e1s1);
		sig2 = RealCalcVec::unpack_hi(e0s0, e1s1);
	#endif

#elif VCP_VEC_TYPE==VCP_VEC_AVX
	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		static const __m256i memoryMask_first_second = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 1<<31);
		const RealCalcVec e0s0 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[0], memoryMask_first_second);
		const RealCalcVec e1s1 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[1], memoryMask_first_second);
		const RealCalcVec e2s2 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[2], memoryMask_first_second);
		const RealCalcVec e3s3 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[3], memoryMask_first_second);
		const RealCalcVec e4s4 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[4], memoryMask_first_second);
		const RealCalcVec e5s5 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[5], memoryMask_first_second);
		const RealCalcVec e6s6 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[6], memoryMask_first_second);
		const RealCalcVec e7s7 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[7], memoryMask_first_second);

		const RealCalcVec e0e1s0s1 = RealCalcVec::unpack_lo(e0s0, e1s1);
		const RealCalcVec e2e3s2s3 = RealCalcVec::unpack_lo(e2s2, e3s3);
		const RealCalcVec e4e5s4s5 = RealCalcVec::unpack_lo(e4s4, e5s5);
		const RealCalcVec e6e7s6s7 = RealCalcVec::unpack_lo(e6s6, e7s7);

		const RealCalcVec e0e1e2e3 = _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(e0e1s0s1), _mm256_castps_pd(e2e3s2s3)));
		const RealCalcVec s0s1s2s3 = _mm256_castpd_ps(_mm256_unpackhi_pd(_mm256_castps_pd(e0e1s0s1), _mm256_castps_pd(e2e3s2s3)));

		const RealCalcVec e4e5e6e7 = _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(e4e5s4s5), _mm256_castps_pd(e6e7s6s7)));
		const RealCalcVec s4s5s6s7 = _mm256_castpd_ps(_mm256_unpackhi_pd(_mm256_castps_pd(e4e5s4s5), _mm256_castps_pd(e6e7s6s7)));

		eps_24 = _mm256_permute2f128_ps(e0e1e2e3, e4e5e6e7, 1<<5);
		sig2 = _mm256_permute2f128_ps(s0s1s2s3, s4s5s6s7, 1<<5);

	#else /* VCP_DPDP */
		static const __m256i memoryMask_first_second = _mm256_set_epi32(0, 0, 0, 0, 1<<31, 0, 1<<31, 0);
		const RealCalcVec e0s0 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[0], memoryMask_first_second);
		const RealCalcVec e1s1 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[1], memoryMask_first_second);
		const RealCalcVec e2s2 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[2], memoryMask_first_second);
		const RealCalcVec e3s3 = RealCalcVec::aligned_load_mask(eps_sigI + 2 * id_j_shifted[3], memoryMask_first_second);

		const RealCalcVec e0e1 = RealCalcVec::unpack_lo(e0s0, e1s1);
		const RealCalcVec s0s1 = RealCalcVec::unpack_hi(e0s0, e1s1);
		const RealCalcVec e2e3 = RealCalcVec::unpack_lo(e2s2, e3s3);
		const RealCalcVec s2s3 = RealCalcVec::unpack_hi(e2s2, e3s3);

		eps_24 = _mm256_permute2f128_pd(e0e1, e2e3, 1<<5);
		sig2 = _mm256_permute2f128_pd(s0s1, s2s3, 1<<5);
	#endif

#elif VCP_VEC_TYPE==VCP_VEC_AVX2//avx
	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		__m256i indices = _mm256_maskload_epi32((const int*)(id_j_shifted), MaskCalcVec::ones());
		indices = _mm256_add_epi32(indices, indices); // only every second...
		eps_24 = _mm256_i32gather_ps(eps_sigI, indices, 4);
		sig2 = _mm256_i32gather_ps(eps_sigI+1, indices, 4);

	#else /* VCP_DPDP */
		__m256i indices = _mm256_maskload_epi64((const long long *)(id_j_shifted), MaskCalcVec::ones());
		indices = _mm256_add_epi64(indices, indices); // only every second...
		eps_24 = _mm256_i64gather_pd(eps_sigI, indices, 8);
		sig2 = _mm256_i64gather_pd(eps_sigI+1, indices, 8);

	#endif /* VCP_PREC */

#elif VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_AVX512F
	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		__m512i indices = _mm512_load_epi32(id_j_shifted);
		indices = _mm512_add_epi32(indices, indices);//only every second...
		eps_24 = _mm512_i32gather_ps(indices, eps_sigI, 4);//eps_sigI+2*id_j[0],eps_sigI+2*id_j[1],...
		sig2 = _mm512_i32gather_ps(indices, eps_sigI+1, 4);//eps_sigI+1+2*id_j[0],eps_sigI+1+2*id_j[1],...
	#else /*VCP_DPDP */
		__m512i indices = _mm512_load_epi64(id_j_shifted);
		indices = _mm512_add_epi64(indices, indices);//only every second...
		eps_24 = _mm512_i64gather_pd(indices, eps_sigI, 8);//eps_sigI+2*id_j[0],eps_sigI+2*id_j[1],...
		sig2 = _mm512_i64gather_pd(indices, eps_sigI+1, 8);//eps_sigI+1+2*id_j[0],eps_sigI+1+2*id_j[1],...
	#endif


#elif VCP_VEC_TYPE==VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER

	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		__m512i indices = _mm512_i32gather_epi32(lookupORforceMask, (const int *) id_j, 4);
		indices = _mm512_add_epi32(indices, indices);//only every second...
		eps_24 = _mm512_i32gather_ps(indices, eps_sigI, 4);//eps_sigI+2*id_j[0],eps_sigI+2*id_j[1],...
		sig2 = _mm512_i32gather_ps(indices, eps_sigI+1, 4);//eps_sigI+1+2*id_j[0],eps_sigI+1+2*id_j[1],...
	#else /*VCP_DPDP*/

		__m256i lookupORforceMask_256i = _mm512_castsi512_si256 (lookupORforceMask);
		__m512i indices = _mm512_i32gather_epi64(lookupORforceMask_256i, (const long long *) id_j, 8);//gather id_j using the indices

		indices = _mm512_add_epi64(indices, indices);//only every second...
		eps_24 = _mm512_i64gather_pd(indices, eps_sigI, 8);//eps_sigI+2*id_j[0],eps_sigI+2*id_j[1],...
		sig2 = _mm512_i64gather_pd(indices, eps_sigI+1, 8);//eps_sigI+1+2*id_j[0],eps_sigI+1+2*id_j[1],...
	#endif /*VCP_PREC*/
#endif
}


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
/**
 * unpacks shift6 and saves it in an DoubleVec (for mic+avx2: use gather)
 * @param shift6 place to store it
 * @param shift6I initial array
 * @param id_j array of displacements
 * @param offset offset of the id_j array
 */
template <class MaskGatherChooser>
static vcp_inline
void unpackShift6(RealCalcVec& shift6, const AlignedArray<vcp_real_calc>& shift6I,
				  const vcp_center_id_t* id_j, const vcp_center_id_t& offset, const vcp_lookupOrMask_vec& lookupORforceMask) {
#if VCP_VEC_TYPE != VCP_VEC_KNL_GATHER
	const vcp_center_id_t* id_j_shifted = id_j + offset;//this is the pointer, to where the stuff is stored.
#endif

#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	shift6 = shift6I[id_j_shifted[0]];

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		const RealCalcVec sh0 = _mm_load_ss(shift6I + id_j_shifted[0]);
		const RealCalcVec sh1 = _mm_load_ss(shift6I + id_j_shifted[1]);
		const RealCalcVec sh2 = _mm_load_ss(shift6I + id_j_shifted[2]);
		const RealCalcVec sh3 = _mm_load_ss(shift6I + id_j_shifted[3]);

		const RealCalcVec sh0sh1 = RealCalcVec::unpack_lo(sh0, sh1);
		const RealCalcVec sh2sh3 = RealCalcVec::unpack_lo(sh2, sh3);

		shift6 = _mm_castpd_ps(_mm_unpacklo_pd(_mm_castps_pd(sh0sh1), _mm_castps_pd(sh2sh3)));

	#else /* VCP_DPDP */
		const RealCalcVec sh0 = _mm_load_sd(shift6I + id_j_shifted[0]);
		const RealCalcVec sh1 = _mm_load_sd(shift6I + id_j_shifted[1]);
		shift6 = RealCalcVec::unpack_lo(sh0, sh1);
	#endif


#elif VCP_VEC_TYPE==VCP_VEC_AVX //avx

	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 1<<31);
		const RealCalcVec sh0 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[0], memoryMask_first);
		const RealCalcVec sh1 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[1], memoryMask_first);
		const RealCalcVec sh2 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[2], memoryMask_first);
		const RealCalcVec sh3 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[3], memoryMask_first);
		const RealCalcVec sh4 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[4], memoryMask_first);
		const RealCalcVec sh5 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[5], memoryMask_first);
		const RealCalcVec sh6 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[6], memoryMask_first);
		const RealCalcVec sh7 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[7], memoryMask_first);

		const RealCalcVec sh0sh1 = RealCalcVec::unpack_lo(sh0, sh1);
		const RealCalcVec sh2sh3 = RealCalcVec::unpack_lo(sh2, sh3);
		const RealCalcVec sh4sh5 = RealCalcVec::unpack_lo(sh4, sh5);
		const RealCalcVec sh6sh7 = RealCalcVec::unpack_lo(sh6, sh7);
		const RealCalcVec sh0sh1sh2sh3 = _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(sh0sh1), _mm256_castps_pd(sh2sh3)));
		const RealCalcVec sh4sh5sh6sh7 = _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(sh4sh5), _mm256_castps_pd(sh6sh7)));
		shift6 = _mm256_permute2f128_ps(sh0sh1sh2sh3, sh4sh5sh6sh7, 1<<5);

	#else /* VCP_PREC == VCP_DPDP */
		static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
		const RealCalcVec sh0 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[0], memoryMask_first);
		const RealCalcVec sh1 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[1], memoryMask_first);
		const RealCalcVec sh2 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[2], memoryMask_first);
		const RealCalcVec sh3 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[3], memoryMask_first);

		const RealCalcVec sh0sh1 = RealCalcVec::unpack_lo(sh0, sh1);
		const RealCalcVec sh2sh3 = RealCalcVec::unpack_lo(sh2, sh3);

		shift6 = _mm256_permute2f128_pd(sh0sh1, sh2sh3, 1<<5);
	#endif

#elif VCP_VEC_TYPE==VCP_VEC_AVX2 //avx2 knows gather

	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		const __m256i indices = _mm256_maskload_epi32((const int*)(id_j_shifted), MaskCalcVec::ones());
		shift6 = _mm256_i32gather_ps(shift6I, indices, 4);

	#else /* VCP_DPDP */
		__m256i indices = _mm256_maskload_epi64((const long long *)(id_j_shifted), MaskCalcVec::ones());
		shift6 = _mm256_i64gather_pd(shift6I, indices, 8);
	#endif



#elif VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_AVX512F

	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		const __m512i indices = _mm512_load_epi32(id_j_shifted);
		shift6 = _mm512_i32gather_ps(indices, shift6I, 4);

	#else /* VCP_DPDP */
		const __m512i indices = _mm512_load_epi64(id_j_shifted);//load id_j, stored continuously
		shift6 = _mm512_i64gather_pd(indices, shift6I, 8);//gather shift6
	#endif


#elif VCP_VEC_TYPE==VCP_VEC_KNL_GATHER or VCP_VEC_TYPE == VCP_VEC_AVX512F_GATHER

	#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
		__m512i indices = _mm512_i32gather_epi32(lookupORforceMask, (const int *) id_j, 4);
		shift6 = _mm512_i32gather_ps(indices, shift6I, 4);

	#else /* VCP_DPDP */
		__m256i lookupORforceMask_256i = _mm512_castsi512_si256 (lookupORforceMask);
		__m512i indices = _mm512_i32gather_epi64(lookupORforceMask_256i, (const long long *) id_j, 8);//gather id_j using the lookupindices
		shift6 = _mm512_i64gather_pd(indices, shift6I, 8);//gather shift6
	#endif

#endif
}
#pragma GCC diagnostic pop

/**
 * sums up values in a and adds the result to *mem_addr
 */
static vcp_inline
void hSum_Add_Store( vcp_real_accum * const mem_addr, const RealAccumVec & a ) {
	RealAccumVec::horizontal_add_and_store(a, mem_addr);
}

static vcp_inline
void load_hSum_Store_Clear(vcp_real_accum * const compact, vcp_real_accum * const wide) {
	RealAccumVec wide_reg = RealAccumVec::aligned_load(wide);
	hSum_Add_Store(compact, wide_reg);
	wide_reg = RealAccumVec::zero();
	wide_reg.aligned_store(wide);
}

/**
 * loads vector from memory location, adds the value to it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be added
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_add_store(vcp_real_accum * const addr, size_t offset, const RealAccumVec& value, const vcp_lookupOrMask_vec& lookupORforceMask){
	RealAccumVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = sum + value;
	MaskGatherChooser::store(addr, offset, sum, lookupORforceMask);
}

/**
 * loads vector from memory location, subtracts the value from it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be subtracted
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_sub_store(vcp_real_accum * const addr, size_t offset, const RealAccumVec& value, const vcp_lookupOrMask_vec& lookupORforceMask){
	RealAccumVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = sum - value;
	MaskGatherChooser::store(addr, offset, sum, lookupORforceMask);
}

template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_fnmadd_store(vcp_real_accum * const addr, size_t offset, const RealAccumVec& value, const RealAccumVec& mult, const vcp_lookupOrMask_vec& lookupORforceMask) {
	RealAccumVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = RealAccumVec::fnmadd(value, mult, sum);
	MaskGatherChooser::store(addr, offset, sum, lookupORforceMask);
}

/**
 * loads vector from memory location, adds the value to it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be added
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_add_store_masked(vcp_real_accum * const addr, size_t offset, const RealAccumVec& value, const vcp_lookupOrMask_vec& lookupORforceMask, const MaskCalcVec mask){
	RealAccumVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = sum + value;
	MaskGatherChooser::storeMasked(addr, offset, sum, lookupORforceMask, mask);
}

/**
 * loads vector from memory location, subtracts the value from it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be subtracted
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_sub_store_masked(vcp_real_accum * const addr, size_t offset, const RealAccumVec& value, const vcp_lookupOrMask_vec& lookupORforceMask, const MaskCalcVec mask){
	RealAccumVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = sum - value;
	MaskGatherChooser::storeMasked(addr, offset, sum, lookupORforceMask, mask);
}


/**
 * \brief Policy class for single cell force calculation.
 */
template<bool ApplyCutoff>
class SingleCellPolicy_ {
public:

	vcp_inline static size_t InitJ (const size_t i)  // needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{  // i: only calculate j>=i
		// however we do a floor for alignment purposes. ->  we have to mark some of the indices to not be computed (this is handled using the InitJ_Mask)
		return vcp_floor_to_vec_size(i);  // this is i if i is divisible by VCP_VEC_SIZE otherwise the next smaller multiple of VCP_VEC_SIZE
	}

	vcp_inline static size_t InitJ2 (const size_t i __attribute__((unused)))  // needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{
#if VCP_VEC_TYPE!=VCP_VEC_KNL_GATHER and VCP_VEC_TYPE!=VCP_VEC_AVX512F_GATHER
		return InitJ(i);
#else
		return 0;
#endif
	}

	vcp_inline static MaskCalcVec InitJ_Mask(const size_t i) {  // calculations only for i onwards.
		return vcp_simd_getInitMask(i);
	}

	vcp_inline static MaskCalcVec GetForceMask(const RealCalcVec& m_r2, const RealCalcVec& rc2, MaskCalcVec& j_mask) {
		MaskCalcVec result = (m_r2 != RealCalcVec::zero()) and j_mask;
		j_mask = MaskCalcVec::ones();

		if (ApplyCutoff == true) {
			 result = (m_r2 < rc2) and result;
		}

		return result;
	}

	vcp_inline static size_t NumDistanceCalculations(const size_t numSoA1, const size_t /*numSoA2*/) {
		return numSoA1 * (numSoA1 - 1) / 2;
	}
}; /* end of class SingleCellPolicy_ */

/**
 * \brief Policy class for cell pair force calculation.
 */
template<bool ApplyCutoff>
class CellPairPolicy_ {
public:

	vcp_inline static size_t InitJ (const size_t /*i*/)
	{
		return 0;
	}
	vcp_inline static size_t InitJ2 (const size_t i)//needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{
		return InitJ(i);
	}

	vcp_inline static MaskCalcVec GetForceMask (const RealCalcVec& m_r2, const RealCalcVec& rc2, MaskCalcVec& /*j_mask*/)
	{
		// Provide a mask with the same logic as used in
		MaskCalcVec result;
		if (ApplyCutoff == true) {
			result = m_r2 < rc2;
		} else {
			result = MaskCalcVec::ones();
		}

		// Note: add m_r2 != 0.0, when (if) Molecules' Sites get split between cells

		return result;
	}

	vcp_inline static MaskCalcVec InitJ_Mask (const size_t /*i*/)
	{
		return MaskCalcVec::ones();//totally unimportant, since not used...
	}

	vcp_inline static size_t NumDistanceCalculations(const size_t numSoA1, const size_t numSoA2) {
		return numSoA1 * numSoA2;
	}
}; /* end of class CellPairPolicy_ */

/**
 * \brief The dist lookup for a molecule and all centers of a type
 */
template<class ForcePolicy, class MaskGatherChooser>
countertype32
static vcp_inline calcDistLookup (const size_t & i_center_idx, const size_t & soa2_num_centers,
		vcp_lookupOrMask_single* const soa2_center_dist_lookup, const vcp_real_calc* const soa2_m_r_x, const vcp_real_calc* const soa2_m_r_y, const vcp_real_calc* const soa2_m_r_z,
		const RealCalcVec & cutoffRadiusSquareD, size_t end_j, const RealCalcVec m1_r_x, const RealCalcVec m1_r_y, const RealCalcVec m1_r_z) {

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	MaskCalcVec initJ_mask = ForcePolicy :: InitJ_Mask(i_center_idx);

	MaskGatherChooser mgc(soa2_center_dist_lookup, j);

	for (; j < end_j; j += VCP_VEC_SIZE) {
		const RealCalcVec m2_r_x = RealCalcVec::aligned_load(soa2_m_r_x + j);
		const RealCalcVec m2_r_y = RealCalcVec::aligned_load(soa2_m_r_y + j);
		const RealCalcVec m2_r_z = RealCalcVec::aligned_load(soa2_m_r_z + j);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		const RealCalcVec m_r2 = RealCalcVec::scal_prod(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const MaskCalcVec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);

		mgc.storeCalcDistLookup(j, forceMask);

	}
	const MaskCalcVec remainderMask = vcp_simd_getRemainderMask(soa2_num_centers);
	if (remainderMask.movemask()) {
		const RealCalcVec m2_r_x = RealCalcVec::aligned_load_mask(soa2_m_r_x + j, remainderMask);
		const RealCalcVec m2_r_y = RealCalcVec::aligned_load_mask(soa2_m_r_y + j, remainderMask);
		const RealCalcVec m2_r_z = RealCalcVec::aligned_load_mask(soa2_m_r_z + j, remainderMask);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		const RealCalcVec m_r2 = RealCalcVec::scal_prod(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const MaskCalcVec forceMask = remainderMask and ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);//AND remainderMask -> set unimportant ones to zero.
		mgc.storeCalcDistLookup(j, forceMask);
	}

	return mgc.getCount();	//do not compute stuff if nothing needs to be computed.

}


#endif /* SIMD_VECTORIZEDCELLPROCESSORHELPERS_H */
