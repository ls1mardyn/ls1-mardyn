#include "SIMD_TYPES.h"

class MaskingChooser {
public:
	inline static size_t getEndloop(const size_t& long_loop, const countertype32& number_calculate /* number of interactions, that are calculated*/) {
		return long_loop;
	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMask(vcp_lookupOrMask_single* dist_lookup, size_t offset){
		return vcp_simd_load(dist_lookup + offset/VCP_INDICES_PER_LOOKUP_SINGLE);
	}

	inline static vcp_double_vec load(const double* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& /*lookup*/) {
		return vcp_simd_load(src + offset);
	}

	inline static void store(double* const addr, const size_t& offset,
			vcp_double_vec& value, const vcp_lookupOrMask_vec& /*lookup*/) {
		vcp_simd_store(addr + offset, value);
	}

	inline static bool computeLoop(const vcp_mask_vec& forceMask) {
		return vcp_simd_movemask(forceMask);
	}

	inline static bool hasRemainder() {
		return false;
	}

	inline static vcp_mask_vec getForceMask(const vcp_mask_vec& forceMask) {
		return forceMask;
	}

};
#if VCP_VEC_TYPE==VCP_VEC_MIC_GATHER
class GatherChooser { //MIC ONLY!!!
public:
	inline static size_t getEndloop(const size_t& long_loop, const countertype32& number_calculate /* number of interactions, that are calculated*/) {
		return vcp_floor_to_vec_size(number_calculate);
	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMask(const vcp_lookupOrMask_single* const dist_lookup, const size_t& offset){
		return _mm512_mask_loadunpackhi_epi32(
			_mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), static_cast<__mmask16>(0x00FF), (dist_lookup + offset)),
			static_cast<__mmask16>(0x00FF),
			dist_lookup + offset + (64 / sizeof(vcp_lookupOrMask_single))
		);
	}

	inline static vcp_lookupOrMask_vec loadLookupOrForceMaskRemainder(const vcp_lookupOrMask_single* const dist_lookup, const size_t& offset, const __mmask8& remainderMask){
		return _mm512_mask_loadunpackhi_epi32(
			_mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), static_cast<__mmask16>(0x00FF) & remainderMask, (dist_lookup + offset)),
			static_cast<__mmask16>(0x00FF) & remainderMask,
			dist_lookup + offset + (64 / sizeof(vcp_lookupOrMask_single))
		);
	}

	inline static vcp_double_vec load(const double* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& lookup) {
		return _mm512_i32logather_pd(lookup, src, 8);
	}

	inline static void store(double* const addr, const size_t& offset,
			vcp_double_vec& value, const vcp_lookupOrMask_vec& lookup) {
		_mm512_i32loscatter_pd(addr, lookup, value, 8);
	}

	inline static void storeMasked(double* const addr, const size_t& offset,
			vcp_double_vec& value, const vcp_lookupOrMask_vec& lookup, const vcp_mask_vec mask) {
		_mm512_mask_i32loscatter_pd(addr, mask, lookup, value, 8);
	}

	inline static bool computeLoop(const vcp_lookupOrMask_vec& forceMask) {
		return true;
	}

	inline static bool hasRemainder() {
		return true;
	}

	inline static vcp_mask_vec getForceMask(const vcp_lookupOrMask_vec& /*forceMask*/) {
		return VCP_SIMD_ONESVM;//compute everything
	}

	inline static __mmask8 getRemainder(const size_t& numberComputations) {
		static const __mmask8 possibleRemainderMasks[8] = { 0x00, 0x01, 0x03, 0x07,
					0x0F, 0x1F, 0x3F, 0x7F };
		return possibleRemainderMasks[numberComputations
				& static_cast<size_t>(7)];
	}
};
#endif

#if VCP_VEC_TYPE!=VCP_VEC_MIC_GATHER
	typedef MaskingChooser MaskGatherC;
#else
	typedef GatherChooser MaskGatherC;
#endif
