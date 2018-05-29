#include "SIMD_TYPES.h"

class MaskingChooser {
private:
	MaskVec compute_molecule=MaskVec::zero();
	vcp_mask_single* const storeCalcDistLookupLocation;
public:
	MaskingChooser(vcp_mask_single* const soa2_center_dist_lookup, size_t /*j*/):
		storeCalcDistLookupLocation(soa2_center_dist_lookup){
	}

	inline int getCount(){
		return compute_molecule.movemask() ? 1 : 0;
	}

	inline void storeCalcDistLookup(size_t j, MaskVec forceMask){
		forceMask.aligned_store(storeCalcDistLookupLocation + j/VCP_INDICES_PER_LOOKUP_SINGLE);
		compute_molecule = compute_molecule or forceMask;
	}

	inline static size_t getEndloop(const size_t& long_loop, const countertype32& /*number_calculate*/ /* number of interactions, that are calculated*/) {
		return long_loop;
	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMask(vcp_lookupOrMask_single* dist_lookup, size_t offset){
		return MaskVec::aligned_load(dist_lookup + offset/VCP_INDICES_PER_LOOKUP_SINGLE);
	}

	inline static RealCalcVec load(const vcp_real_calc* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& /*lookup*/) {
		return RealCalcVec::aligned_load(src + offset);
	}

	inline static void store(vcp_real_calc* const addr, const size_t& offset,
			RealCalcVec& value, const vcp_lookupOrMask_vec& /*lookup*/) {
		value.aligned_store(addr + offset);
	}

	inline static bool computeLoop(const MaskVec& forceMask) {
		return forceMask.movemask();
	}

	inline static bool hasRemainder() {
		return false;
	}

	inline static MaskVec getForceMask(const MaskVec& forceMask) {
		return forceMask;
	}

};
#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER or VCP_VEC_TYPE==VCP_VEC_KNL_GATHER
class GatherChooser { //MIC ONLY!!!
private:
	__m512i indices;
	vcp_lookupOrMask_single* const storeCalcDistLookupLocation;
	int counter = 0;
public:
	GatherChooser(vcp_lookupOrMask_single* const soa2_center_dist_lookup, size_t j):
		storeCalcDistLookupLocation(soa2_center_dist_lookup)
	{
		static const __m512i first_indices = _mm512_set_epi32(
			0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
			0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00
		);
		indices = _mm512_mask_add_epi32(_mm512_setzero_epi32(), static_cast<__mmask16>(0x00FF), first_indices, _mm512_set1_epi32(j));
	}

	inline void storeCalcDistLookup(size_t j, MaskVec forceMask){
		#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
			_mm512_mask_packstorelo_epi32(storeCalcDistLookupLocation + counter, static_cast<__mmask16>(forceMask), indices);//these two lines are an unaligned store
			_mm512_mask_packstorehi_epi32(storeCalcDistLookupLocation + counter + (64 / sizeof(vcp_lookupOrMask_single)), static_cast<__mmask16>(forceMask), indices);//these two lines are an unaligned store
		#else
			_mm512_mask_compressstoreu_epi32(storeCalcDistLookupLocation + counter, static_cast<__mmask16>(forceMask), indices);
		#endif

		static const __m512i eight = _mm512_set_epi32(
			0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
			0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08, 0x08
		);

		indices = _mm512_add_epi32(indices, eight);
		counter += _popcnt64(forceMask);
	}
	inline int getCount(){
		return counter;
	}
	inline static size_t getEndloop(const size_t& long_loop, const countertype32& number_calculate /* number of interactions, that are calculated*/) {
		return vcp_floor_to_vec_size(number_calculate);
	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMask(const vcp_lookupOrMask_single* const dist_lookup, const size_t& offset){
		#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
			return _mm512_mask_loadunpackhi_epi32(
				_mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), static_cast<__mmask16>(0x00FF), (dist_lookup + offset)),
				static_cast<__mmask16>(0x00FF),
				dist_lookup + offset + (64 / sizeof(vcp_lookupOrMask_single))
			);
		#else
			return _mm512_maskz_loadu_epi32(static_cast<__mmask16>(0x00FF), dist_lookup + offset);
		#endif

	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMaskRemainder(const vcp_lookupOrMask_single* const dist_lookup, const size_t& offset, const __mmask8& remainderMask){
		#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
			return _mm512_mask_loadunpackhi_epi32(
				_mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), static_cast<__mmask16>(0x00FF) & remainderMask, (dist_lookup + offset)),
				static_cast<__mmask16>(0x00FF) & remainderMask,
				dist_lookup + offset + (64 / sizeof(vcp_lookupOrMask_single))
			);
		#else
			return _mm512_maskz_loadu_epi32(static_cast<__mmask16>(0x00FF) & remainderMask, dist_lookup + offset);
		#endif
	}

	inline static RealCalcVec load(const vcp_real_calc* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& lookup) {
		#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
			return _mm512_i32logather_pd(lookup, src, 8);
		#else
			__m256i lookup_256i = _mm512_castsi512_si256 (lookup);
			return _mm512_i32gather_pd(lookup_256i, src, 8);
		#endif
	}

	inline static void store(vcp_real_calc* const addr, const size_t& offset,
			RealCalcVec& value, const vcp_lookupOrMask_vec& lookup) {
		#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
			_mm512_i32loscatter_pd(addr, lookup, value, 8);
		#else
			__m256i lookup_256i = _mm512_castsi512_si256 (lookup);
			_mm512_i32scatter_pd(addr, lookup_256i, value, 8);
		#endif
	}

	inline static void storeMasked(vcp_real_calc* const addr, const size_t& offset,
			RealCalcVec& value, const vcp_lookupOrMask_vec& lookup, const MaskVec mask) {
		#if VCP_VEC_TYPE==VCP_VEC_KNC_GATHER
			_mm512_mask_i32loscatter_pd(addr, mask, lookup, value, 8);
		#else
			__m256i lookup_256i = _mm512_castsi512_si256 (lookup);
			_mm512_mask_i32scatter_pd(addr, mask, lookup_256i, value, 8);
		#endif
	}

	inline static bool computeLoop(const vcp_lookupOrMask_vec& forceMask) {
		return true;
	}

	inline static bool hasRemainder() {
		return true;
	}

	inline static MaskVec getForceMask(const vcp_lookupOrMask_vec& /*forceMask*/) {
		return MaskVec::ones();//compute everything
	}

	inline static __mmask8 getRemainder(const size_t& numberComputations) {
		static const __mmask8 possibleRemainderMasks[8] = { 0x00, 0x01, 0x03, 0x07,
					0x0F, 0x1F, 0x3F, 0x7F };
		return possibleRemainderMasks[numberComputations
				& static_cast<size_t>(7)];
	}
};
#endif

#if VCP_VEC_TYPE!=VCP_VEC_KNC_GATHER and VCP_VEC_TYPE!=VCP_VEC_KNL_GATHER
	typedef MaskingChooser MaskGatherC;
#else
	typedef GatherChooser MaskGatherC;
#endif
