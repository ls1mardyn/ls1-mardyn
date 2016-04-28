/*
 * NonBlockingMPIMultiStepHandler.cpp
 *
 *  Created on: Apr 28, 2016
 *      Author: seckler
 */

#include "NonBlockingMPIMultiStepHandler.h"
#include <assert.h>

using Log::global_log;
using namespace std;

NonBlockingMPIMultiStepHandler::NonBlockingMPIMultiStepHandler(
		Timer* decompositionTimer, Timer* computationTimer,
		DomainDecompMPIBase* domainDecomposition,
		ParticleContainer* moleculeContainer, Domain* domain,
		CellProcessor* cellProcessor) :
		NonBlockingMPIHandlerBase(decompositionTimer, computationTimer,
				domainDecomposition, moleculeContainer, domain, cellProcessor) {
	// just calling the constructor of the base class is enough

}

NonBlockingMPIMultiStepHandler::~NonBlockingMPIMultiStepHandler() {
	// nothing to be done
}

void NonBlockingMPIMultiStepHandler::performComputation() {
	int stageCount = _domainDecomposition->getNonBlockingStageCount();

	assert(stageCount > 0);

	for (unsigned int i = 0; i < stageCount; ++i) {
		_decompositionTimer->start();
		_domainDecomposition->prepareNonBlockingStage(forceRebalancing,
				_moleculeContainer, _domain, i);
		_decompositionTimer->stop();

		_computationTimer->start();
		// Force calculation and other pair interaction related computations
		global_log->debug() << "Traversing pairs" << std::endl;
		_moleculeContainer->traversePartialInnermostCells(*_cellProcessor, i,
				stageCount);
		_computationTimer->stop();

		_decompositionTimer->start();
		_domainDecomposition->finishNonBlockingStage(forceRebalancing,
				_moleculeContainer, _domain, i);
		_decompositionTimer->stop();
	}
	_decompositionTimer->start();
	_moleculeContainer->updateBoundaryAndHaloMoleculeCaches();//update the caches of the other molecules (non-inner cells)
	_decompositionTimer->stop();

	_computationTimer->start();
	// Force calculation and other pair interaction related computations
	global_log->debug() << "Traversing pairs" << std::endl;
	_moleculeContainer->traverseNonInnermostCells(*_cellProcessor);
	_computationTimer->stop();
}

void NonBlockingMPIMultiStepHandler::initBalanceAndExchange(
		bool forceRebalancing) {

	_decompositionTimer->start();
	_domainDecomposition->balanceAndExchangeInitNonBlocking(forceRebalancing,
			_moleculeContainer, _domain);

	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateInnerMoleculeCaches(); // only the caches of the innermost molecules have to be updated.

	_decompositionTimer->stop();
}

