/*
 * P2MCellProcessor.h
 *
 *  Created on: Feb 10, 2015
 *      Author: tchipev
 */

#ifndef P2MCELLPROCESSOR_H_
#define P2MCELLPROCESSOR_H_

#include "particleContainer/adapter/CellProcessor.h"
#include "utils/Timer.h"

namespace bhfmm {

class PseudoParticleContainer;

class P2MCellProcessor: public CellProcessor {
public:
	P2MCellProcessor(PseudoParticleContainer * pseudoParticleContainer);
	~P2MCellProcessor();

	virtual double processSingleMolecule(Molecule* m1, ParticleCell& cell2) {return 0.0;}

	void initTraversal(const size_t numCells) {}
	void preprocessCell(ParticleCell& cell) {}
	void processCellPair(ParticleCell& cell1, ParticleCell& cell2) {}
	void processCell(ParticleCell& cell);
	void postprocessCell(ParticleCell& cell) {}
	void endTraversal() {}

	void printTimers();

private:
	PseudoParticleContainer* const _pseudoParticleContainer;
	Timer _P2MTimer;
};

} /* namespace bhfmm */

#endif /* P2MCELLPROCESSOR_H_ */
