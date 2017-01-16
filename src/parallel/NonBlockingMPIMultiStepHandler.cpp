/*
 * NonBlockingMPIMultiStepHandler.cpp
 *
 *  Created on: Apr 28, 2016
 *      Author: seckler
 */

#include "NonBlockingMPIMultiStepHandler.h"
#include <assert.h>
#include "utils/Logger.h"
#include "utils/Timer.h"
#include "DomainDecompMPIBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/CellProcessor.h"

using Log::global_log;
using namespace std;

NonBlockingMPIMultiStepHandler::NonBlockingMPIMultiStepHandler(Timer* decompositionTimer,
		Timer* computationTimer, Timer* forceCalculationTimer,
		DomainDecompMPIBase* domainDecomposition,
		ParticleContainer* moleculeContainer, Domain* domain,
		CellProcessor* cellProcessor) :
				NonBlockingMPIHandlerBase(decompositionTimer, computationTimer, forceCalculationTimer,
						domainDecomposition, moleculeContainer, domain, cellProcessor) {
	// just calling the constructor of the base class is enough
}

NonBlockingMPIMultiStepHandler::~NonBlockingMPIMultiStepHandler() {
	// nothing to be done
}

void NonBlockingMPIMultiStepHandler::performComputation() {
	int stageCount = _domainDecomposition->getNonBlockingStageCount();

	assert(stageCount > 0);
	_forceCalculationTimer->start();
	_cellProcessor->initTraversal();
	_forceCalculationTimer->stop();
	for (unsigned int i = 0; i < static_cast<unsigned int>(stageCount); ++i) {
		_decompositionTimer->start();
		_domainDecomposition->prepareNonBlockingStage(false,
				_moleculeContainer, _domain, i);
		_decompositionTimer->stop();

		_computationTimer->start();
		// Force calculation and other pair interaction related computations
		global_log->debug() << "Traversing innermost cells" << std::endl;
		_forceCalculationTimer->start();
		_moleculeContainer->traversePartialInnermostCells(*_cellProcessor, i,
				stageCount);
		_forceCalculationTimer->stop();
		_computationTimer->stop();

		_decompositionTimer->start();
		_domainDecomposition->finishNonBlockingStage(false,
				_moleculeContainer, _domain, i);
		_decompositionTimer->stop();
	}

	_decompositionTimer->start();
	_moleculeContainer->updateBoundaryAndHaloMoleculeCaches();  // update the caches of the other molecules (non-inner cells)
	_decompositionTimer->stop();

	_computationTimer->start();
	// remaining force calculation and other pair interaction related computations
	global_log->debug() << "Traversing non-innermost cells" << std::endl;
	_forceCalculationTimer->start();
	_moleculeContainer->traverseNonInnermostCells(*_cellProcessor);
	_cellProcessor->endTraversal();
	_forceCalculationTimer->stop();
	_computationTimer->stop();
}

void NonBlockingMPIMultiStepHandler::initBalanceAndExchange(
		bool forceRebalancing) {

	assert(!forceRebalancing);

	_decompositionTimer->start();
	_domainDecomposition->balanceAndExchangeInitNonBlocking(forceRebalancing,
			_moleculeContainer, _domain);

	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateInnerMoleculeCaches(); // only the caches of the innermost molecules have to be updated.

	_decompositionTimer->stop();
}

