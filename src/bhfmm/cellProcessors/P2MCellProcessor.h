/*
 * P2MCellProcessor.h
 *
 *  Created on: Feb 10, 2015
 *      Author: tchipev
 */

#ifndef P2MCELLPROCESSOR_H_
#define P2MCELLPROCESSOR_H_

#include "bhfmm/cellProcessors/SimpleCellProcessor.h"
#include "utils/Timer.h"
#include <stdlib.h>

namespace bhfmm {

class PseudoParticleContainer;

class P2MCellProcessor: public SimpleCellProcessor {
public:
	P2MCellProcessor(PseudoParticleContainer * pseudoParticleContainer);
	~P2MCellProcessor();

	void initTraversal() {}
	void processCell(ParticleCellPointers& cell);
	void endTraversal() {}

	void printTimers();

private:
	PseudoParticleContainer* const _pseudoParticleContainer;
	Timer _P2MTimer;
};

} /* namespace bhfmm */

#endif /* P2MCELLPROCESSOR_H_ */
