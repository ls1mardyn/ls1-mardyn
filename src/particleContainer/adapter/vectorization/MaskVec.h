/*
 * MaskVec.h
 *
 *  Created on: 14 Nov 2016
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_MASKVEC_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_MASKVEC_H_

#include "./SIMD_TYPES.h"

namespace vcp {

template<typename FloatOrDouble>
class MaskVec {
};

} /* namespace vcp */

// SPECIALIZATIONS: suggested to view via vimdiff
#include "MaskVecFloat.h"
#include "MaskVecDouble.h"



#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_MASKVEC_H_ */
