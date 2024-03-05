/*
 * L2PCellProcessor.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: tchipev
 */

#include "L2PCellProcessor.h"
#include "Simulation.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "bhfmm/containers/PseudoParticleContainer.h"


namespace bhfmm {

L2PCellProcessor::L2PCellProcessor(
		PseudoParticleContainer * pseudoParticleContainer) :
		_pseudoParticleContainer(pseudoParticleContainer) {
#ifdef ENABLE_MPI
	global_simulation->timers()->setOutputString("L2P_CELL_PROCESSOR_L2P", "FMM: Time spent in L2P ");
	//global_simulation->setSyncTimer("L2P_CELL_PROCESSOR_L2P", false); //it is per default false
#endif
}

L2PCellProcessor::~L2PCellProcessor() {
}

void L2PCellProcessor::initTraversal() {
//	using std::cout;
//	using std::endl;
//	Domain* domain = global_simulation->getDomain();
//	std::cout << "L2P init: LocalUpot     " << domain->getLocalUpot() << std::endl;
//	std::cout << "L2P init: LocalVirial   " << domain->getLocalVirial() << std::endl;
//	std::cout << "L2P init: LocalP_xx     " << domain->getLocalP_xx() << std::endl;
//	std::cout << "L2P init: LocalP_yy     " << domain->getLocalP_yy() << std::endl;
//	std::cout << "L2P init: LocalP_zz     " << domain->getLocalP_zz() << std::endl;
	global_simulation->timers()->start("L2P_CELL_PROCESSOR_L2P");
}

void L2PCellProcessor::processCell(ParticleCellPointers& cell) {
	if (!cell.isHaloCell()) {
		_pseudoParticleContainer->processFarField(cell);
	}
}

void L2PCellProcessor::printTimers() {
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	int numprocs = domainDecomp.getNumProcs();
	int myrank = domainDecomp.getRank();
	for (int i = 0; i < numprocs; i++) {
		if (i == myrank) {
			std::cout << "rank: " << myrank << std::endl;
			std::cout << "\t\t" << global_simulation->timers()->getTime("L2P_CELL_PROCESSOR_L2P") << "\t\t" << "s in L2P" << std::endl;
			global_simulation->timers()->print("L2P_CELL_PROCESSOR_L2P");
		}
		domainDecomp.barrier();
	}
}

void L2PCellProcessor::endTraversal() {
	global_simulation->timers()->stop("L2P_CELL_PROCESSOR_L2P");
}

} /* namespace bhfmm */
