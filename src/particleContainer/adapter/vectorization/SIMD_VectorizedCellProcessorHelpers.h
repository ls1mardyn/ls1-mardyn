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

#if VCP_VEC_TYPE==VCP_VEC_AVX2 or \
	VCP_VEC_TYPE==VCP_VEC_KNC or \
	VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or \
	VCP_VEC_TYPE==VCP_VEC_KNL or \
	VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
static_assert (sizeof(size_t) == 8, "Code assumes, that sizeof(size_t) is 8. Contact SCCS developers if this fails.");
#endif /* sizeof(size_t) */

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
void unpackEps24Sig2(RealCalcVec& eps_24, RealCalcVec& sig2, const DoubleArray& eps_sigI,
		const size_t* const id_j, const size_t& offset, const vcp_lookupOrMask_vec& lookupORforceMask __attribute__((unused))) {

#if VCP_VEC_TYPE != VCP_VEC_KNC_GATHER and VCP_VEC_TYPE != VCP_VEC_KNL_GATHER
	const size_t* id_j_shifted = id_j + offset;//this is the pointer, to where the stuff is stored.
#endif

#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	eps_24 = eps_sigI[2 * id_j_shifted[0]];
	sig2 = eps_sigI[2 * id_j_shifted[0] + 1];

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	const RealCalcVec e0s0 = RealCalcVec::aligned_load(eps_sigI + 2 * id_j_shifted[0]);
	const RealCalcVec e1s1 = RealCalcVec::aligned_load(eps_sigI + 2 * id_j_shifted[1]);
	eps_24 = RealCalcVec::unpack_lo(e0s0, e1s1);
	sig2 = RealCalcVec::unpack_hi(e0s0, e1s1);

#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2//avx
	//TODO: add gather for AVX2 (here only)
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

#elif VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNL
	__m512i indices = _mm512_load_epi64(id_j_shifted);
	indices = _mm512_add_epi64(indices, indices);//only every second...
	eps_24 = _mm512_i64gather_pd(indices, eps_sigI, 8);//eps_sigI+2*id_j[0],eps_sigI+2*id_j[1],...
	sig2 = _mm512_i64gather_pd(indices, eps_sigI+1, 8);//eps_sigI+1+2*id_j[0],eps_sigI+1+2*id_j[1],...

#elif VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER

	#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
		__m512i indices = _mm512_i32logather_epi64(lookupORforceMask, id_j, 8);//gather id_j using the indices
	#else
		__m256i lookupORforceMask_256i = _mm512_castsi512_si256 (lookupORforceMask);
		__m512i indices = _mm512_i32gather_epi64(lookupORforceMask_256i, id_j, 8);//gather id_j using the indices
	#endif

	indices = _mm512_add_epi64(indices, indices);//only every second...
	eps_24 = _mm512_i64gather_pd(indices, eps_sigI, 8);//eps_sigI+2*id_j[0],eps_sigI+2*id_j[1],...
	sig2 = _mm512_i64gather_pd(indices, eps_sigI+1, 8);//eps_sigI+1+2*id_j[0],eps_sigI+1+2*id_j[1],...
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
void unpackShift6(RealCalcVec& shift6, const DoubleArray& shift6I,
		const size_t* id_j, const size_t& offset, const vcp_lookupOrMask_vec& lookupORforceMask) {
#if VCP_VEC_TYPE != VCP_VEC_KNC_GATHER and VCP_VEC_TYPE != VCP_VEC_KNL_GATHER
	const size_t* id_j_shifted = id_j + offset;//this is the pointer, to where the stuff is stored.
#endif

#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	shift6 = shift6I[id_j_shifted[0]];

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	const RealCalcVec sh1 = _mm_load_sd(shift6I + id_j_shifted[0]);
	const RealCalcVec sh2 = _mm_load_sd(shift6I + id_j_shifted[1]);
	shift6 = RealCalcVec::unpack_lo(sh1, sh2);

#elif VCP_VEC_TYPE==VCP_VEC_AVX //avx
	static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
	const RealCalcVec sh0 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[0], memoryMask_first);
	const RealCalcVec sh1 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[1], memoryMask_first);
	const RealCalcVec sh2 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[2], memoryMask_first);
	const RealCalcVec sh3 = RealCalcVec::aligned_load_mask(shift6I + id_j_shifted[3], memoryMask_first);

	const RealCalcVec sh0sh1 = RealCalcVec::unpack_lo(sh0, sh1);
	const RealCalcVec sh2sh3 = RealCalcVec::unpack_lo(sh2, sh3);

	shift6 = _mm256_permute2f128_pd(sh0sh1, sh2sh3, 1<<5);

#elif VCP_VEC_TYPE==VCP_VEC_AVX2 //avx2 knows gather
	const __m256i indices = _mm256_maskload_epi64((const long long*)(id_j_shifted), MaskVec::ones());
	shift6 = _mm256_i64gather_pd(shift6I, indices, 8);

#elif VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNL
	const __m512i indices = _mm512_load_epi64(id_j_shifted);//load id_j, stored continuously
	shift6 = _mm512_i64gather_pd(indices, shift6I, 8);//gather shift6

#elif VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER

	#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
		__m512i indices = _mm512_i32logather_epi64(lookupORforceMask, id_j, 8);//gather id_j using the lookupindices
	#else
		__m256i lookupORforceMask_256i = _mm512_castsi512_si256 (lookupORforceMask);
		__m512i indices = _mm512_i32gather_epi64(lookupORforceMask_256i, id_j, 8);//gather id_j using the lookupindices
	#endif

	shift6 = _mm512_i64gather_pd(indices, shift6I, 8);//gather shift6
#endif
}
#pragma GCC diagnostic pop


#if VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
static vcp_inline
double horizontal_add_256 (__m256d a) {
    __m256d t1 = _mm256_hadd_pd(a,a);
    __m128d t2 = _mm256_extractf128_pd(t1,1);
    __m128d t3 = _mm_add_sd(_mm256_castpd256_pd128(t1),t2);
    return _mm_cvtsd_f64(t3);
}
#endif

/**
 * sums up values in a and adds the result to *mem_addr
 */
static vcp_inline
void hSum_Add_Store( double * const mem_addr, const RealCalcVec & a ) {
#if VCP_VEC_TYPE==VCP_NOVEC
	(*mem_addr) += a; //there is just one value of a, so no second sum needed.
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	_mm_store_sd(
		mem_addr,
		_mm_add_sd(
			_mm_hadd_pd(a, a),
			_mm_load_sd(mem_addr)
		)
	);
#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2
	static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
	const RealCalcVec a_t1 = _mm256_permute2f128_pd(a, a, 0x1);
	const RealCalcVec a_t2 = _mm256_hadd_pd(a, a_t1);
	const RealCalcVec a_t3 = _mm256_hadd_pd(a_t2, a_t2);
	_mm256_maskstore_pd(
		mem_addr,
		memoryMask_first,
			a_t3 + RealCalcVec::aligned_load_mask(mem_addr, memoryMask_first)
	);
#elif VCP_VEC_TYPE==VCP_VEC_KNC or VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
	*mem_addr += _mm512_reduce_add_pd(a);

#elif VCP_VEC_TYPE==VCP_VEC_KNL or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
	// NOTE: separate, because only the Intel compiler provides _mm512_reduce_add_pd
	__m256d low = _mm512_castpd512_pd256(a);
	__m256d high = _mm512_extractf64x4_pd(a, 1);
	*mem_addr += horizontal_add_256(low + high);
#endif
}

static vcp_inline
void load_hSum_Store_Clear(double * const compact, double * const wide) {
	RealCalcVec wide_reg = RealCalcVec::aligned_load(wide);
	hSum_Add_Store(compact, wide_reg);
	wide_reg = RealCalcVec::zero();
	wide_reg.aligned_store(wide);
}

/**
 * loads vector from memory location, adds the value to it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be added
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_add_store(double * const addr, size_t offset, const RealCalcVec& value, const vcp_lookupOrMask_vec& lookupORforceMask){
	RealCalcVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
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
void vcp_simd_load_sub_store(double * const addr, size_t offset, const RealCalcVec& value, const vcp_lookupOrMask_vec& lookupORforceMask){
	RealCalcVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = sum - value;
	MaskGatherChooser::store(addr, offset, sum, lookupORforceMask);
}

/**
 * loads vector from memory location, adds the value to it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be added
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_add_store_masked(double * const addr, size_t offset, const RealCalcVec& value, const vcp_lookupOrMask_vec& lookupORforceMask, const MaskVec mask){
	RealCalcVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
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
void vcp_simd_load_sub_store_masked(double * const addr, size_t offset, const RealCalcVec& value, const vcp_lookupOrMask_vec& lookupORforceMask, const MaskVec mask){
	RealCalcVec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = sum - value;
	MaskGatherChooser::storeMasked(addr, offset, sum, lookupORforceMask, mask);
}


/**
 * \brief Policy class for single cell force calculation.
 */
template<bool ApplyCutoff>
class SingleCellPolicy_ {
public:

	inline static size_t InitJ (const size_t i)  // needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{  // i: only calculate j>=i
		// however we do a floor for alignment purposes. ->  we have to mark some of the indices to not be computed (this is handled using the InitJ_Mask)
		return vcp_floor_to_vec_size(i);  // this is i if i is divisible by VCP_VEC_SIZE otherwise the next smaller multiple of VCP_VEC_SIZE
	}

	inline static size_t InitJ2 (const size_t i __attribute__((unused)))  // needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{
#if VCP_VEC_TYPE!=VCP_VEC_KNC_GATHER and VCP_VEC_TYPE!=VCP_VEC_KNL_GATHER
		return InitJ(i);
#else
		return 0;
#endif
	}

	inline static MaskVec InitJ_Mask(const size_t i) {  // calculations only for i onwards.
		return vcp_simd_getInitMask(i);
	}

	inline static MaskVec GetForceMask(const RealCalcVec& m_r2, const RealCalcVec& rc2, MaskVec& j_mask) {
		MaskVec result = (m_r2 != RealCalcVec::zero()) and j_mask;
		j_mask = MaskVec::ones();

		if (ApplyCutoff == true) {
			 result = (m_r2 < rc2) and result;
		}

		return result;
	}
}; /* end of class SingleCellPolicy_ */

/**
 * \brief Policy class for cell pair force calculation.
 */
template<bool ApplyCutoff>
class CellPairPolicy_ {
public:

	inline static size_t InitJ (const size_t /*i*/)
	{
		return 0;
	}
	inline static size_t InitJ2 (const size_t i)//needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{
		return InitJ(i);
	}

	inline static MaskVec GetForceMask (const RealCalcVec& m_r2, const RealCalcVec& rc2, MaskVec& /*j_mask*/)
	{
		// Provide a mask with the same logic as used in
		MaskVec result;
		if (ApplyCutoff == true) {
			result = m_r2 < rc2;
		} else {
			result = MaskVec::ones();
		}

		// Note: add m_r2 != 0.0, when (if) Molecules' Sites get split between cells

		return result;
	}

	inline static MaskVec InitJ_Mask (const size_t /*i*/)
	{
		return MaskVec::ones();//totally unimportant, since not used...
	}
}; /* end of class CellPairPolicy_ */

/**
 * \brief The dist lookup for a molecule and all centers of a type
 */
template<class ForcePolicy, class MaskGatherChooser>
countertype32
static inline calcDistLookup (const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
		vcp_lookupOrMask_single* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z,
		const RealCalcVec & cutoffRadiusSquareD, size_t end_j, const RealCalcVec m1_r_x, const RealCalcVec m1_r_y, const RealCalcVec m1_r_z) {

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	MaskVec initJ_mask = ForcePolicy :: InitJ_Mask(i_center_idx);

	MaskGatherChooser mgc(soa2_center_dist_lookup, j);

	for (; j < end_j; j += VCP_VEC_SIZE) {
		const RealCalcVec m2_r_x = RealCalcVec::aligned_load(soa2_m_r_x + j);
		const RealCalcVec m2_r_y = RealCalcVec::aligned_load(soa2_m_r_y + j);
		const RealCalcVec m2_r_z = RealCalcVec::aligned_load(soa2_m_r_z + j);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		const RealCalcVec m_r2 = RealCalcVec::scal_prod(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const MaskVec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);

		mgc.storeCalcDistLookup(j, forceMask);

	}
	const MaskVec remainderMask = vcp_simd_getRemainderMask(soa2_num_centers);
	if (remainderMask.movemask()) {
		const RealCalcVec m2_r_x = RealCalcVec::aligned_load_mask(soa2_m_r_x + j, remainderMask);
		const RealCalcVec m2_r_y = RealCalcVec::aligned_load_mask(soa2_m_r_y + j, remainderMask);
		const RealCalcVec m2_r_z = RealCalcVec::aligned_load_mask(soa2_m_r_z + j, remainderMask);

		const RealCalcVec m_dx = m1_r_x - m2_r_x;
		const RealCalcVec m_dy = m1_r_y - m2_r_y;
		const RealCalcVec m_dz = m1_r_z - m2_r_z;

		const RealCalcVec m_r2 = RealCalcVec::scal_prod(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const MaskVec forceMask = remainderMask and ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);//AND remainderMask -> set unimportant ones to zero.
		mgc.storeCalcDistLookup(j, forceMask);
	}

	return mgc.getCount();	//do not compute stuff if nothing needs to be computed.

}


#endif /* SIMD_VECTORIZEDCELLPROCESSORHELPERS_H */
