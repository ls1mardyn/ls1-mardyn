/*
 * L2PCellProcessor.h
 *
 *  Created on: Feb 10, 2015
 *      Author: tchipev
 */

#ifndef L2PCELLPROCESSOR_H_
#define L2PCELLPROCESSOR_H_
#include "bhfmm/cellProcessors/SimpleCellProcessor.h"
#include <stdlib.h>

namespace bhfmm {

class PseudoParticleContainer;

class L2PCellProcessor: public SimpleCellProcessor {
public:
	L2PCellProcessor(PseudoParticleContainer * pseudoParticleContainer);
	~L2PCellProcessor();

	void initTraversal();
	void processCell(ParticleCellPointers& cell);
	void endTraversal();

	void printTimers();

private:
	PseudoParticleContainer* const _pseudoParticleContainer;
};

} /* namespace bhfmm */

#endif /* L2PCELLPROCESSOR_H_ */
