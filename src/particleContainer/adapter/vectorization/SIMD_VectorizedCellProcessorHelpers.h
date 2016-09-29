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
void unpackEps24Sig2(vcp_double_vec& eps_24, vcp_double_vec& sig2, const DoubleArray& eps_sigI,
		const size_t* const id_j, const size_t& offset, const vcp_lookupOrMask_vec& lookupORforceMask __attribute__((unused))){
#if VCP_VEC_TYPE != VCP_VEC_MIC_GATHER
	const size_t* id_j_shifted = id_j + offset;//this is the pointer, to where the stuff is stored.
#endif
#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	eps_24 = eps_sigI[2 * id_j_shifted[0]];
	sig2 = eps_sigI[2 * id_j_shifted[0] + 1];

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	const vcp_double_vec e0s0 = vcp_simd_load(eps_sigI + 2 * id_j_shifted[0]);
	const vcp_double_vec e1s1 = vcp_simd_load(eps_sigI + 2 * id_j_shifted[1]);
	eps_24 = vcp_simd_unpacklo(e0s0, e1s1);
	sig2 = vcp_simd_unpackhi(e0s0, e1s1);

#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2//avx
	static const __m256i memoryMask_first_second = _mm256_set_epi32(0, 0, 0, 0, 1<<31, 0, 1<<31, 0);
	const vcp_double_vec e0s0 = vcp_simd_maskload(eps_sigI + 2 * id_j_shifted[0], memoryMask_first_second);
	const vcp_double_vec e1s1 = vcp_simd_maskload(eps_sigI + 2 * id_j_shifted[1], memoryMask_first_second);
	const vcp_double_vec e2s2 = vcp_simd_maskload(eps_sigI + 2 * id_j_shifted[2], memoryMask_first_second);
	const vcp_double_vec e3s3 = vcp_simd_maskload(eps_sigI + 2 * id_j_shifted[3], memoryMask_first_second);

	const vcp_double_vec e0e1 = vcp_simd_unpacklo(e0s0, e1s1);
	const vcp_double_vec s0s1 = vcp_simd_unpackhi(e0s0, e1s1);
	const vcp_double_vec e2e3 = vcp_simd_unpacklo(e2s2, e3s3);
	const vcp_double_vec s2s3 = vcp_simd_unpackhi(e2s2, e3s3);

	eps_24 = _mm256_permute2f128_pd(e0e1, e2e3, 1<<5);
	sig2 = _mm256_permute2f128_pd(s0s1, s2s3, 1<<5);
#elif VCP_VEC_TYPE==VCP_VEC_MIC //mic knows gather, too
		#define BUILD_BUG_ON1212(condition) ((void)sizeof(char[1 - 2*!!(condition)]))
		BUILD_BUG_ON1212(sizeof(size_t) % 8);//check whether size_t is of size 8...
	__m512i indices = _mm512_load_epi64(id_j_shifted);
	indices = _mm512_add_epi64(indices, indices);//only every second...
	eps_24 = _mm512_i64gather_pd(indices, eps_sigI, 8);//eps_sigI+2*id_j[0],eps_sigI+2*id_j[1],...
	sig2 = _mm512_i64gather_pd(indices, eps_sigI+1, 8);//eps_sigI+1+2*id_j[0],eps_sigI+1+2*id_j[1],...
#elif VCP_VEC_TYPE==VCP_VEC_MIC_GATHER
	__m512i indices = _mm512_i32logather_epi64(lookupORforceMask, id_j, 8);//gather id_j using the indices
	indices = _mm512_add_epi64(indices, indices);//only every second...
	eps_24 = _mm512_i64gather_pd(indices, eps_sigI, 8);//eps_sigI+2*id_j[0],eps_sigI+2*id_j[1],...
	sig2 = _mm512_i64gather_pd(indices, eps_sigI+1, 8);//eps_sigI+1+2*id_j[0],eps_sigI+1+2*id_j[1],...
#endif
}


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
/**
 * unpacks shift6 and saves it in an vcp_double_vec (for mic+avx2: use gather)
 * @param shift6 place to store it
 * @param shift6I initial array
 * @param id_j array of displacements
 * @param offset offset of the id_j array
 */
template <class MaskGatherChooser>
static vcp_inline
void unpackShift6(vcp_double_vec& shift6, const DoubleArray& shift6I,
		const size_t* id_j, const size_t& offset, const vcp_lookupOrMask_vec& lookupORforceMask){
#if VCP_VEC_TYPE != VCP_VEC_MIC_GATHER
	const size_t* id_j_shifted = id_j + offset;//this is the pointer, to where the stuff is stored.
#endif
#if VCP_VEC_TYPE==VCP_NOVEC //novec comes first. For NOVEC no specific types are specified -- use build in ones.
	shift6 = shift6I[id_j_shifted[0]];

#elif VCP_VEC_TYPE==VCP_VEC_SSE3 //sse3
	const vcp_double_vec sh1 = _mm_load_sd(shift6I + id_j_shifted[0]);
	const vcp_double_vec sh2 = _mm_load_sd(shift6I + id_j_shifted[1]);
	shift6 = vcp_simd_unpacklo(sh1, sh2);

#elif VCP_VEC_TYPE==VCP_VEC_AVX //avx
	static const __m256i memoryMask_first = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 1<<31, 0);
	const vcp_double_vec sh0 = vcp_simd_maskload(shift6I + id_j_shifted[0], memoryMask_first);
	const vcp_double_vec sh1 = vcp_simd_maskload(shift6I + id_j_shifted[1], memoryMask_first);
	const vcp_double_vec sh2 = vcp_simd_maskload(shift6I + id_j_shifted[2], memoryMask_first);
	const vcp_double_vec sh3 = vcp_simd_maskload(shift6I + id_j_shifted[3], memoryMask_first);

	const vcp_double_vec sh0sh1 = vcp_simd_unpacklo(sh0, sh1);
	const vcp_double_vec sh2sh3 = vcp_simd_unpacklo(sh2, sh3);

	shift6 = _mm256_permute2f128_pd(sh0sh1, sh2sh3, 1<<5);
#elif VCP_VEC_TYPE==VCP_VEC_AVX2 //avx2 knows gather
		#define BUILD_BUG_ON1212(condition) ((void)sizeof(char[1 - 2*!!(condition)]))
		BUILD_BUG_ON1212(sizeof(size_t) % 8);//check whether size_t is of size 8...
	const __m256i indices = _mm256_maskload_epi64((const long long*)(id_j_shifted), VCP_SIMD_ONESVM);
	shift6 = _mm256_i64gather_pd(shift6I, indices, 8);
#elif VCP_VEC_TYPE==VCP_VEC_MIC //mic knows gather, too
		#define BUILD_BUG_ON1212(condition) ((void)sizeof(char[1 - 2*!!(condition)]))
		BUILD_BUG_ON1212(sizeof(size_t) % 8);//check whether size_t is of size 8...
	const __m512i indices = _mm512_load_epi64(id_j_shifted);//load id_j, stored continuously
	shift6 = _mm512_i64gather_pd(indices, shift6I, 8);//gather shift6
#elif VCP_VEC_TYPE==VCP_VEC_MIC_GATHER
	__m512i indices = _mm512_i32logather_epi64(lookupORforceMask, id_j, 8);//gather id_j using the lookupindices
	shift6 = _mm512_i64gather_pd(indices, shift6I, 8);//gather shift6
#endif
}
#pragma GCC diagnostic pop


/**
 * sums up values in a and adds the result to *mem_addr
 */
static vcp_inline
void hSum_Add_Store( double * const mem_addr, const vcp_double_vec & a ) {
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
	const vcp_double_vec a_t1 = _mm256_permute2f128_pd(a, a, 0x1);
	const vcp_double_vec a_t2 = _mm256_hadd_pd(a, a_t1);
	const vcp_double_vec a_t3 = _mm256_hadd_pd(a_t2, a_t2);
	_mm256_maskstore_pd(
		mem_addr,
		memoryMask_first,
		vcp_simd_add(
			a_t3,
			vcp_simd_maskload(mem_addr, memoryMask_first)
		)
	);
#elif VCP_VEC_TYPE==VCP_VEC_MIC or VCP_VEC_TYPE==VCP_VEC_MIC_GATHER
	*mem_addr += _mm512_reduce_add_pd(a);
#endif
}

static vcp_inline
void load_hSum_Store_Clear(double * const compact, double * const wide) {
	vcp_double_vec wide_reg = vcp_simd_load(wide);
	hSum_Add_Store(compact, wide_reg);
	wide_reg = vcp_simd_zerov();
	vcp_simd_store(wide, wide_reg);
}

/**
 * loads vector from memory location, adds the value to it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be added
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_add_store(double * const addr, size_t offset, const vcp_double_vec& value, const vcp_lookupOrMask_vec& lookupORforceMask){
	vcp_double_vec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = vcp_simd_add(sum, value);
	MaskGatherChooser::store(addr, offset, sum, lookupORforceMask);
}

/**
 * loads vector from memory location, subtracts the value from it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be subtracted
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_sub_store(double * const addr, size_t offset, const vcp_double_vec& value, const vcp_lookupOrMask_vec& lookupORforceMask){
	vcp_double_vec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = vcp_simd_sub(sum, value);
	MaskGatherChooser::store(addr, offset, sum, lookupORforceMask);
}

/**
 * loads vector from memory location, adds the value to it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be added
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_add_store_masked(double * const addr, size_t offset, const vcp_double_vec& value, const vcp_lookupOrMask_vec& lookupORforceMask, const vcp_mask_vec mask){
	vcp_double_vec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = vcp_simd_add(sum, value);
	MaskGatherChooser::storeMasked(addr, offset, sum, lookupORforceMask, mask);
}

/**
 * loads vector from memory location, subtracts the value from it and saves the combined result.
 * @param addr memory address where value should be loaded from and stored to
 * @param value value that should be subtracted
 */
template <class MaskGatherChooser>
static vcp_inline
void vcp_simd_load_sub_store_masked(double * const addr, size_t offset, const vcp_double_vec& value, const vcp_lookupOrMask_vec& lookupORforceMask, const vcp_mask_vec mask){
	vcp_double_vec sum = MaskGatherChooser::load(addr, offset, lookupORforceMask);
	sum = vcp_simd_sub(sum, value);
	MaskGatherChooser::storeMasked(addr, offset, sum, lookupORforceMask, mask);
}


/**
 * \brief Policy class for single cell force calculation.
 */
class SingleCellPolicy_ {
public:

	inline static size_t InitJ (const size_t i)  // needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{  // i: only calculate j>=i
		// however we do a floor for alignment purposes. ->  we have to mark some of the indices to not be computed (this is handled using the InitJ_Mask)
		return vcp_floor_to_vec_size(i);  // this is i if i is divisible by VCP_VEC_SIZE otherwise the next smaller multiple of VCP_VEC_SIZE
	}

	inline static size_t InitJ2 (const size_t i __attribute__((unused)))  // needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{
#if VCP_VEC_TYPE!=VCP_VEC_MIC_GATHER
		return InitJ(i);
#else
		return 0;
#endif
	}

	inline static vcp_mask_vec InitJ_Mask(const size_t i) {  // calculations only for i onwards.
		return vcp_simd_getInitMask(i);
	}

	inline static vcp_mask_vec GetForceMask(const vcp_double_vec& m_r2, const vcp_double_vec& rc2,
			vcp_mask_vec& j_mask) {
		vcp_mask_vec result = vcp_simd_and(vcp_simd_and(vcp_simd_lt(m_r2, rc2), vcp_simd_neq(m_r2, VCP_SIMD_ZEROV)),
				j_mask);
		j_mask = VCP_SIMD_ONESVM;
		return result;
	}



}; /* end of class SingleCellPolicy_ */

/**
 * \brief Policy class for cell pair force calculation.
 */
class CellPairPolicy_ {
public:

	inline static bool Condition(double m_r2, double rc2)
	{
		// Because we have 2 different cells, no 2 centers on the same
		// molecule can form a pair.
		return m_r2 < rc2;
	}

	inline static size_t InitJ (const size_t /*i*/)
	{
		return 0;
	}
	inline static size_t InitJ2 (const size_t i)//needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
	{
		return InitJ(i);
	}

	inline static vcp_mask_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2, vcp_mask_vec& /*j_mask*/)
	{
		// Provide a mask with the same logic as used in
		// bool Condition(double m_r2, double rc2)
		return vcp_simd_lt(m_r2, rc2);
	}

	inline static vcp_mask_vec InitJ_Mask (const size_t /*i*/)
	{
		return VCP_SIMD_ONESVM;//totally unimportant, since not used...
	}

}; /* end of class CellPairPolicy_ */

/**
 * \brief The dist lookup for a molecule and all centers of a type
 * \author Robert Hajda
 */
template<class ForcePolicy, class MaskGatherChooser>
countertype32
static inline calcDistLookup (const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
		vcp_lookupOrMask_single* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z,
		const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z) {

	size_t j = ForcePolicy :: InitJ(i_center_idx);
	vcp_mask_vec initJ_mask = ForcePolicy :: InitJ_Mask(i_center_idx);

	MaskGatherChooser mgc(soa2_center_dist_lookup, j);

	for (; j < end_j; j += VCP_VEC_SIZE) {
		const vcp_double_vec m2_r_x = vcp_simd_load(soa2_m_r_x + j);
		const vcp_double_vec m2_r_y = vcp_simd_load(soa2_m_r_y + j);
		const vcp_double_vec m2_r_z = vcp_simd_load(soa2_m_r_z + j);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = vcp_simd_scalProd(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const vcp_mask_vec forceMask = ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask);

		mgc.storeCalcDistLookup(j, forceMask);

	}
	const vcp_mask_vec remainderMask = vcp_simd_getRemainderMask(soa2_num_centers);
	if (vcp_simd_movemask(remainderMask)) {
		const vcp_double_vec m2_r_x = vcp_simd_maskload(soa2_m_r_x + j, remainderMask);
		const vcp_double_vec m2_r_y = vcp_simd_maskload(soa2_m_r_y + j, remainderMask);
		const vcp_double_vec m2_r_z = vcp_simd_maskload(soa2_m_r_z + j, remainderMask);

		const vcp_double_vec m_dx = m1_r_x - m2_r_x;
		const vcp_double_vec m_dy = m1_r_y - m2_r_y;
		const vcp_double_vec m_dz = m1_r_z - m2_r_z;

		const vcp_double_vec m_r2 = vcp_simd_scalProd(m_dx, m_dy, m_dz, m_dx, m_dy, m_dz);

		const vcp_mask_vec forceMask = vcp_simd_and(remainderMask, ForcePolicy::GetForceMask(m_r2, cutoffRadiusSquareD, initJ_mask));//AND remainderMask -> set unimportant ones to zero.
		mgc.storeCalcDistLookup(j, forceMask);
	}

	return mgc.getCount();	//do not compute stuff if nothing needs to be computed.

}


#endif /* SIMD_VECTORIZEDCELLPROCESSORHELPERS_H */
