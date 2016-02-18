/*
 * FastMultipoleMethod.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: tchipev
 */

#include "FastMultipoleMethod.h"
#include "Simulation.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "bhfmm/containers/UniformPseudoParticleContainer.h"
#include "bhfmm/containers/AdaptivePseudoParticleContainer.h"

using Log::global_log;
using std::endl;

namespace bhfmm {

FastMultipoleMethod::FastMultipoleMethod(double globalDomainLength[3],
		double bBoxMin[3], double bBoxMax[3], double LJCellLength[3],
		unsigned LJSubdivisionFactor, int orderOfExpansions, bool periodic, bool adaptive) :
		_order(orderOfExpansions), _wellSeparated(1), _pseudoParticleContainer(0),
		_P2PProcessor(0), _P2MProcessor(0), _L2PProcessor(0) {

	if (LJSubdivisionFactor != 1 and
		LJSubdivisionFactor != 2 and
		LJSubdivisionFactor != 4 and
		LJSubdivisionFactor != 8) {
		global_log->error() << "Fast Multipole Method: bad subdivision factor:" << LJSubdivisionFactor << endl;
		global_log->error() << "expected 1,2,4 or 8" << endl;
		exit(5);
	}
	global_log->info() << "Fast Multipole Method: each LJ cell will be subdivided in "
			<< pow(LJSubdivisionFactor, 3) << " cells for electrostatic calculations in FMM" << endl;

	_P2PProcessor = new VectorizedChargeP2PCellProcessor(*(global_simulation->getDomain()));

	if (not adaptive) {
		_pseudoParticleContainer = new UniformPseudoParticleContainer(
				globalDomainLength, bBoxMin, bBoxMax, LJCellLength,
				LJSubdivisionFactor, _order, periodic);

	} else {
		// TODO: Debugging in Progress!
		//int threshold = 100;
		//_pseudoParticleContainer = new AdaptivePseudoParticleContainer(
		//		globalDomainLength, threshold, _order, periodic);
		_pseudoParticleContainer = new AdaptivePseudoParticleContainer(
				globalDomainLength, _order, LJCellLength, LJSubdivisionFactor,
				periodic);
	}

	_P2MProcessor = new P2MCellProcessor(_pseudoParticleContainer);
	_L2PProcessor = new L2PCellProcessor(_pseudoParticleContainer);
  
}

FastMultipoleMethod::~FastMultipoleMethod() {
	delete _pseudoParticleContainer;
	delete _P2PProcessor;
	delete _P2MProcessor;
	delete _L2PProcessor;
}

void FastMultipoleMethod::computeElectrostatics(ParticleContainer* ljContainer) {
	// build
	_pseudoParticleContainer->build(ljContainer);

	// clear expansions
	_pseudoParticleContainer->clear();

	// P2M, M2P
	_pseudoParticleContainer->upwardPass(_P2MProcessor);

	// M2L, P2P
	_pseudoParticleContainer->horizontalPass(_P2PProcessor);

	// L2L, L2P
	_pseudoParticleContainer->downwardPass(_L2PProcessor);

}

void FastMultipoleMethod::printTimers() {
	_P2PProcessor->printTimers();
	_P2MProcessor->printTimers();
	_L2PProcessor->printTimers();
	_pseudoParticleContainer->printTimers();
}

} /* namespace bhfmm */

