#pragma once

#include "SIMD_TYPES.h"

class MaskingChooser {
private:
	MaskCalcVec compute_molecule=MaskCalcVec::zero();
	vcp_mask_single* const storeCalcDistLookupLocation;
public:
	MaskingChooser(vcp_mask_single* const soa2_center_dist_lookup, size_t /*j*/):
		storeCalcDistLookupLocation(soa2_center_dist_lookup){
	}

	inline int getCount(){
		return compute_molecule.movemask() ? 1 : 0;
	}

	inline void storeCalcDistLookup(size_t j, MaskCalcVec forceMask){
		forceMask.aligned_store(storeCalcDistLookupLocation + j/VCP_INDICES_PER_LOOKUP_SINGLE);
		compute_molecule = compute_molecule or forceMask;
	}

	inline static size_t getEndloop(const size_t& long_loop, const countertype32& /*number_calculate*/ /* number of interactions, that are calculated*/) {
		return long_loop;
	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMask(vcp_lookupOrMask_single* dist_lookup, size_t offset){
		return MaskCalcVec::aligned_load(dist_lookup + offset/VCP_INDICES_PER_LOOKUP_SINGLE);
	}

	inline static RealCalcVec load(const vcp_real_calc* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& /*lookup*/) {
		return RealCalcVec::aligned_load(src + offset);
	}

	inline static void store(vcp_real_calc* const addr, const size_t& offset,
			RealCalcVec& value, const vcp_lookupOrMask_vec& /*lookup*/) {
		value.aligned_store(addr + offset);
	}
#if VCP_PREC == VCP_SPDP
	inline static RealAccumVec load(const vcp_real_accum* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& /*lookup*/) {
		return RealAccumVec::aligned_load(src + offset);
	}

	inline static void store(vcp_real_accum* const addr, const size_t& offset,
			RealAccumVec& value, const vcp_lookupOrMask_vec& /*lookup*/) {
		value.aligned_store(addr + offset);
	}
#endif

	inline static bool computeLoop(const MaskCalcVec& forceMask) {
		return forceMask.movemask();
	}

	inline static bool hasRemainder() {
		return false;
	}

	inline static MaskCalcVec getForceMask(const MaskCalcVec& forceMask) {
		return forceMask;
	}

};
#if VCP_VEC_TYPE==VCP_VEC_KNL_GATHER or VCP_VEC_TYPE==VCP_VEC_AVX512F_GATHER
class GatherChooser { //scatter needed, i.e. AVX512
private:
	__m512i indices;
	vcp_lookupOrMask_single* const storeCalcDistLookupLocation;
	int counter = 0;
public:
	GatherChooser(vcp_lookupOrMask_single* const soa2_center_dist_lookup, size_t j):
		storeCalcDistLookupLocation(soa2_center_dist_lookup)
	{
		static const __m512i first_indices = _mm512_set_epi32(
			0x0f, 0x0e, 0x0d, 0x0c, 0x0b, 0x0a, 0x09, 0x08,
			0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00
		);

		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			indices = _mm512_add_epi32(first_indices, _mm512_set1_epi32(j));
		#else /* VCP_DPDP */
			indices = _mm512_mask_add_epi32(_mm512_setzero_epi32(), static_cast<__mmask16>(0x00FF), first_indices, _mm512_set1_epi32(j));
		#endif
	}

	inline void storeCalcDistLookup(size_t j, MaskCalcVec forceMask){
		_mm512_mask_compressstoreu_epi32(storeCalcDistLookupLocation + counter, static_cast<__mmask16>(forceMask), indices);

		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			static const __m512i advance = _mm512_set_epi32(
				0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10,
				0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10
			);
		#else /* VCP_DPDP */
			static const __m512i advance = _mm512_set_epi32(
				0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
				0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08
			);
		#endif

		//advance the indices to the next 8 (DP) or 16 (SP) values
		indices = _mm512_add_epi32(indices, advance);
		counter += __builtin_popcount(forceMask);
	}
	inline int getCount(){
		return counter;
	}
	inline static size_t getEndloop(const size_t& long_loop, const countertype32& number_calculate /* number of interactions, that are calculated*/) {
		return vcp_floor_to_vec_size(number_calculate);
	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMask(const vcp_lookupOrMask_single* const dist_lookup, const size_t& offset){
		return _mm512_maskz_loadu_epi32(static_cast<__mmask16>(0xFFFF), dist_lookup + offset);
	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMaskRemainder(const vcp_lookupOrMask_single* const dist_lookup, const size_t& offset, const __mmask16& remainderMask){
		return _mm512_maskz_loadu_epi32(static_cast<__mmask16>(0xFFFF) & remainderMask, dist_lookup + offset);
	}

	inline static RealCalcVec load(const vcp_real_calc* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& lookup) {
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			return _mm512_i32gather_ps(lookup, src, 4);
		#else
			__m256i lookup_256i = _mm512_castsi512_si256 (lookup);
			return _mm512_i32gather_pd(lookup_256i, src, 8);
		#endif
	}

	inline static void store(vcp_real_calc* const addr, const size_t& offset,
			RealCalcVec& value, const vcp_lookupOrMask_vec& lookup) {
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			_mm512_i32scatter_ps(addr, lookup, value, 4);
		#else
			__m256i lookup_256i = _mm512_castsi512_si256 (lookup);
			_mm512_i32scatter_pd(addr, lookup_256i, value, 8);
		#endif
	}

	inline static void storeMasked(vcp_real_calc* const addr, const size_t& offset,
			RealCalcVec& value, const vcp_lookupOrMask_vec& lookup, const MaskCalcVec mask) {
		#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
			_mm512_mask_i32scatter_ps(addr, mask, lookup, value, 4);
		#else
			__m256i lookup_256i = _mm512_castsi512_si256 (lookup);
			_mm512_mask_i32scatter_pd(addr, mask, lookup_256i, value, 8);
		#endif
	}


#if VCP_PREC == VCP_SPDP
	inline static RealAccumVec load(const vcp_real_accum* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& lookup) {
		return RealAccumVec::gather_load(src, offset, lookup);
	}

	inline static void store(vcp_real_accum* const addr, const size_t& offset,
			RealAccumVec& value, const vcp_lookupOrMask_vec& lookup) {
		value.scatter_store(addr, offset, lookup);
	}

	inline static void storeMasked(vcp_real_accum* const addr, const size_t& offset,
			RealAccumVec& value, const vcp_lookupOrMask_vec& lookup, const MaskCalcVec mask) {
		value.scatter_store_mask(addr, offset, lookup, mask);
	}
#endif


	inline static bool computeLoop(const vcp_lookupOrMask_vec& forceMask) {
		return true;
	}

	inline static bool hasRemainder() {
		return true;
	}

	inline static MaskCalcVec getForceMask(const vcp_lookupOrMask_vec& /*forceMask*/) {
		return MaskCalcVec::ones();//compute everything
	}

#if VCP_PREC == VCP_SPSP or VCP_PREC == VCP_SPDP
	inline static __mmask16 getRemainder(const size_t& numberComputations) {
		static const __mmask16 possibleRemainderMasks[16] = { 0x0000, 0x0001, 0x0003, 0x0007, 0x000F, 0x001F, 0x003F, 0x007F,
															  0x00FF, 0x01FF, 0x03FF, 0x07FF, 0x0FFF, 0x1FFF, 0x3FFF, 0x7FFF };
		return possibleRemainderMasks[numberComputations & static_cast<size_t>(15)];
	}
#else /* VCP_DPDP */
	inline static __mmask8 getRemainder(const size_t& numberComputations) {
		static const __mmask8 possibleRemainderMasks[8] = { 0x00, 0x01, 0x03, 0x07,
					0x0F, 0x1F, 0x3F, 0x7F };
		return possibleRemainderMasks[numberComputations
				& static_cast<size_t>(7)];
	}
#endif
};
#endif

#if VCP_VEC_TYPE!=VCP_VEC_KNL_GATHER and VCP_VEC_TYPE!=VCP_VEC_AVX512F_GATHER
	typedef MaskingChooser MaskGatherC;
#else
	typedef GatherChooser MaskGatherC;
#endif

class CountUnmasked_MGC {
private:
	countertype32 _numUnmasked;
public:
	vcp_inline CountUnmasked_MGC(vcp_lookupOrMask_single* const /*soa2_center_dist_lookup*/, size_t /*j*/)
	: _numUnmasked(0) {}

	vcp_inline int getCount() {
		return _numUnmasked;
	}

	vcp_inline void storeCalcDistLookup(size_t j, MaskCalcVec forceMask) {
		_numUnmasked += forceMask.countUnmasked();
	}
};
