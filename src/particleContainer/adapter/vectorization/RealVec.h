/*
 * RealVec.h
 *
 *  Created on: 10 Nov 2016
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVEC_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVEC_H_

#include "SIMD_TYPES.h"
#include "MaskVec.h"

#include <fstream>
#include <cmath>

#include "utils/Logger.h"

namespace vcp {

template<typename FloatOrDouble>
class RealVec {
}; /* class RealVec */

} /* namespace vcp */

// SPECIALIZATIONS: suggested to view via vimdiff
#include "RealVecFloat.h"
#include "RealVecDouble.h"



#endif /* SRC_PARTICLECONTAINER_ADAPTER_VECTORIZATION_REALVEC_H_ */
