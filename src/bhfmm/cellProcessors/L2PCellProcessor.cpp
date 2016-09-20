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
	_L2PTimer.set_sync(false);
#endif
}

L2PCellProcessor::~L2PCellProcessor() {
}

void L2PCellProcessor::initTraversal() {
//	using std::cout;
//	using std::endl;
//	Domain* domain = global_simulation->getDomain();
//	cout << "L2P init: LocalUpot     " << domain->getLocalUpot() << endl;
//	cout << "L2P init: LocalVirial   " << domain->getLocalVirial() << endl;
//	cout << "L2P init: LocalP_xx     " << domain->getLocalP_xx() << endl;
//	cout << "L2P init: LocalP_yy     " << domain->getLocalP_yy() << endl;
//	cout << "L2P init: LocalP_zz     " << domain->getLocalP_zz() << endl;
}

void L2PCellProcessor::processCell(ParticleCell& cell) {
	if (!cell.isHaloCell()) {
		_L2PTimer.start();
		_pseudoParticleContainer->processFarField(cell);
		_L2PTimer.stop();
	}
}

void L2PCellProcessor::printTimers() {
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	int numprocs = domainDecomp.getNumProcs();
	int myrank = domainDecomp.getRank();
	for (int i = 0; i < numprocs; i++) {
		if (i == myrank) {
			std::cout << "rank: " << myrank << std::endl;
			std::cout << "\t\t" << _L2PTimer.get_etime() << "\t\t" << "s in L2P"
					<< std::endl;
		}
		domainDecomp.barrier();
	}

}

void L2PCellProcessor::endTraversal() {
}

} /* namespace bhfmm */
