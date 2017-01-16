/*
 * P2MCellProcessor.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: tchipev
 */

#include "P2MCellProcessor.h"
#include "Simulation.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "bhfmm/containers/PseudoParticleContainer.h"
#include "bhfmm/containers/ParticleCellPointers.h"


namespace bhfmm {

P2MCellProcessor::P2MCellProcessor(
		PseudoParticleContainer * pseudoParticleContainer) :
		_pseudoParticleContainer(pseudoParticleContainer) {
#ifdef ENABLE_MPI
	_P2MTimer.set_sync(false);
#endif
}

P2MCellProcessor::~P2MCellProcessor() {
}

void P2MCellProcessor::processCell(ParticleCellPointers& cell) {
	if (!cell.isHaloCell()) {
		_P2MTimer.start();
		_pseudoParticleContainer->processMultipole(cell);
		_P2MTimer.stop();
	}
}

void P2MCellProcessor::printTimers() {
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	int numprocs = domainDecomp.getNumProcs();
	int myrank = domainDecomp.getRank();
	for (int i = 0; i < numprocs; i++) {
		if (i == myrank) {
			std::cout << "rank: " << myrank << std::endl;
			std::cout << "\t\t" << _P2MTimer.get_etime() << "\t\t" << "s in P2M" << std::endl;
		}
		domainDecomp.barrier();
	}

}

} /* namespace bhfmm */

