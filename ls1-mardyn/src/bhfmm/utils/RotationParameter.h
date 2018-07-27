/*
 * RotationParameter.h
 *
 *  Created on: Jul 7, 2015
 *      Author: uwe
 */

#ifndef ROTATIONPARAMETER_H_
#define ROTATIONPARAMETER_H_

#include "bhfmm/utils/WignerMatrix.h"

namespace bhfmm {

class RotationParams {
public:
//	~RotationParams(){ delete SinCos; }
	WignerMatrix W[2];
	double* SinCos;
};

}  // namespace bhfmm



#endif /* ROTATIONPARAMETER_H_ */
