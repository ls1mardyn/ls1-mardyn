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

	virtual double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) {return 0.0;}
    virtual int countNeighbours(Molecule* /*m1*/, ParticleCell& /*cell2*/, double /*RR*/) { exit(0); return 0; }

	void initTraversal() {}
	void processCell(ParticleCell& cell);
	void endTraversal() {}

	void printTimers();

private:
	PseudoParticleContainer* const _pseudoParticleContainer;
	Timer _P2MTimer;
};

} /* namespace bhfmm */

#endif /* P2MCELLPROCESSOR_H_ */
