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
#include "utils/xmlfileUnits.h"

using Log::global_log;
using std::endl;

namespace bhfmm {

FastMultipoleMethod::~FastMultipoleMethod() {
	delete _pseudoParticleContainer;
	delete _P2PProcessor;
	delete _P2MProcessor;
	delete _L2PProcessor;
}

void FastMultipoleMethod::readXML(XMLfileUnits& xmlconfig) {

	xmlconfig.getNodeValue("orderOfExpansions", _order);
	global_log->info() << "FastMultipoleMethod: orderOfExpansions: " << _order << endl;

	xmlconfig.getNodeValue("LJCellSubdivisionFactor", _LJCellSubdivisionFactor);
	global_log->info() << "FastMultipoleMethod: LJCellSubdivisionFactor: " << _LJCellSubdivisionFactor << endl;

	xmlconfig.getNodeValue("adaptiveContainer", _adaptive);
	if (_adaptive == 1) {
		global_log->warning() << "FastMultipoleMethod: adaptiveContainer is not debugged yet and certainly delivers WRONG results!" << endl;
		global_log->warning() << "Unless you are in the process of debugging this container, please stop the simulation and restart with the uniform one" << endl;
	} else {
		global_log->info() << "FastMultipoleMethod: UniformPseudoParticleSelected " << endl;
	}

	xmlconfig.getNodeValue("systemIsPeriodic", _periodic);
	if (_periodic == 0) {
		global_log->warning() << "FastMultipoleMethod: periodicity is turned off!" << endl;
	} else {
		global_log->info() << "FastMultipoleMethod: Periodicity is on." << endl;
	}
}

void FastMultipoleMethod::setParameters(unsigned LJSubdivisionFactor,
		int orderOfExpansions, bool periodic, bool adaptive) {
	_LJCellSubdivisionFactor = LJSubdivisionFactor;
	_order = orderOfExpansions;
	_periodic = periodic;
	_adaptive = adaptive;
}

void FastMultipoleMethod::init(double globalDomainLength[3], double bBoxMin[3],
		double bBoxMax[3], double LJCellLength[3]) {

	if (_LJCellSubdivisionFactor != 1
			and _LJCellSubdivisionFactor != 2
			and _LJCellSubdivisionFactor != 4
			and _LJCellSubdivisionFactor != 8) {
		global_log->error() << "Fast Multipole Method: bad subdivision factor:"
				<< _LJCellSubdivisionFactor << endl;
		global_log->error() << "expected 1,2,4 or 8" << endl;
		exit(5);
	}
	global_log->info()
			<< "Fast Multipole Method: each LJ cell will be subdivided in "
			<< pow(_LJCellSubdivisionFactor, 3)
			<< " cells for electrostatic calculations in FMM" << endl;

	_P2PProcessor = new VectorizedChargeP2PCellProcessor(
			*(global_simulation->getDomain()));
	if (not _adaptive) {
		_pseudoParticleContainer = new UniformPseudoParticleContainer(
				globalDomainLength, bBoxMin, bBoxMax, LJCellLength,
				_LJCellSubdivisionFactor, _order, _periodic);

	} else {
		// TODO: Debugging in Progress!
		//int threshold = 100;
		_pseudoParticleContainer = new AdaptivePseudoParticleContainer(
				globalDomainLength, _order, LJCellLength,
				_LJCellSubdivisionFactor, _periodic);
	}

	_P2MProcessor = new P2MCellProcessor(_pseudoParticleContainer);
	_L2PProcessor = new L2PCellProcessor(_pseudoParticleContainer);

}

void FastMultipoleMethod::computeElectrostatics(ParticleContainer* ljContainer) {
	// build
	_pseudoParticleContainer->build(ljContainer);

	// clear expansions
	_pseudoParticleContainer->clear();

	// P2M, M2P
	_pseudoParticleContainer->upwardPass(_P2MProcessor);

	// M2L, P2P
	if (_adaptive) {
		_P2PProcessor->initTraversal(2);
	}
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

