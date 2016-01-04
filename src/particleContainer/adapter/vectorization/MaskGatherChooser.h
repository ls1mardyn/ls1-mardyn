#include "SIMD_TYPES.h"

class MaskingChooser {
public:
	inline static size_t getEndloop(const size_t& full_number) {
		return full_number;
	} //TODO: add second parameter for gather

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
private:
	static const __mmask8 possibleRemainderMasks[8] = { 0x00, 0x01, 0x03, 0x07,
			0x0F, 0x1F, 0x3F, 0x7F };
public:
	inline static size_t getEndloop(const size_t& full_number) {
		return full_number;
	} //TODO: add second parameter for gather

	inline static vcp_double_vec load(const double* const src,
			const size_t& offset, const vcp_lookupOrMask_vec& lookup) {
		__m512i indices = _mm512_load_epi64(lookup + offset);
		return _mm512_i64gather_pd(indices, src, 8);
	}

	inline static void store(double* const addr, const size_t& offset,
			vcp_double_vec& value, const vcp_lookupOrMask_vec& lookup) {
		vcp_simd_store(addr + offset, value);
		const __m512i indices = _mm512_load_epi64(lookup + offset);
		_mm512_i64scatter_pd(addr, indices, value, 8);
	}

	inline static bool computeLoop(const vcp_lookupOrMask_vec& forceMask) {
		return true;
	}

	inline static bool hasRemainder() {
		return false;
	}

	inline static vcp_mask_vec getForceMask(const size_t*& /*forceMask*/) {
		return VCP_SIMD_ONESVM;
	}

	inline static __mmask8 getRemainder(size_t numberComputations) {
		return possibleRemainderMasks[numberComputations
				& static_cast<size_t>(7)];
	}
};
#endif
