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
	global_simulation->setOutputString("P2M_CELL_PROCESSOR_P2M", "FMM: Time spent in P2M ");
	//global_simulation->setSyncTimer("P2M_CELL_PROCESSOR_P2M", false); //it is per default false
#endif
}

P2MCellProcessor::~P2MCellProcessor() {
}

void P2MCellProcessor::processCell(ParticleCellPointers& cell) {
	if (!cell.isHaloCell()) {
		global_simulation->startTimer("P2M_CELL_PROCESSOR_P2M");
		_pseudoParticleContainer->processMultipole(cell);
		global_simulation->stopTimer("P2M_CELL_PROCESSOR_P2M");
	}
}

void P2MCellProcessor::printTimers() {
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	int numprocs = domainDecomp.getNumProcs();
	int myrank = domainDecomp.getRank();
	for (int i = 0; i < numprocs; i++) {
		if (i == myrank) {
			std::cout << "rank: " << myrank << std::endl;
			std::cout << "\t\t" << global_simulation->getTime("P2M_CELL_PROCESSOR_P2M") << "\t\t" << "s in P2M" << std::endl;
			//global_simulation->printTimer("P2M_CELL_PROCESSOR_P2M");
		}
		domainDecomp.barrier();
	}
}

} /* namespace bhfmm */

