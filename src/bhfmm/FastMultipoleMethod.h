/*
 * FastMultipoleMethod.h
 *
 *  Created on: Feb 7, 2015
 *      Author: tchipev
 */

#ifndef FASTMULTIPOLEMETHOD_H_
#define FASTMULTIPOLEMETHOD_H_

#include "bhfmm/containers/PseudoParticleContainer.h"
#include "particleContainer/ParticleContainer.h"

namespace bhfmm {

class FastMultipoleMethod {
public:
	FastMultipoleMethod(double globalDomainLength[3], double localBBoxMin[3],
			double localBBoxMax[3], double LJCellLength[3], unsigned LJSubdivisionFactor,
			int orderOfExpansions,
			bool periodic = true, bool adaptive=false); // don't try this at home!

	~FastMultipoleMethod();

	void computeElectrostatics(ParticleContainer * ljContainer);

	void printTimers();

private:
	int _order;
	int _wellSeparated;

	PseudoParticleContainer * _pseudoParticleContainer;

	VectorizedChargeP2PCellProcessor *_P2PProcessor;
	P2MCellProcessor *_P2MProcessor;
	L2PCellProcessor *_L2PProcessor;

};

} /* namespace bhfmm */

#endif /* FASTMULTIPOLEMETHOD_H_ */
