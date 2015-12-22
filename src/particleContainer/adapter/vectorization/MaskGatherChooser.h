#include "SIMD_TYPES.h"

class MaskingChooser{
public:
	inline static size_t getEndloop(const size_t& full_number){return full_number;}//TODO: add second parameter for gather
	inline static vcp_double_vec load(const double* const src, const size_t& offset){return vcp_simd_load(src + offset);}//TODO: add parameter for gather (index array)
	inline static void store(double* const addr, const size_t& offset, vcp_double_vec& value){vcp_simd_store(addr + offset, value);}//TODO: add parameter for gather (index array)
	inline static bool computeLoop(const vcp_mask_vec& forceMask){return vcp_simd_movemask(forceMask);}//TODO: true for gatherchooser

};
