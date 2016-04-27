/*
 * NonBlockingMPIBase.cpp
 *
 *  Created on: Apr 27, 2016
 *      Author: seckler
 */

#include "NonBlockingMPIHandlerBase.h"
#include "utils/Timer.h"
#include "utils/Logger.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "particleContainer/adapter/CellProcessor.h"

using Log::global_log;

NonBlockingMPIHandlerBase::NonBlockingMPIHandlerBase(Timer* decompositionTimer,
		Timer* computationTimer, DomainDecompBase* domainDecomposition,
		ParticleContainer* moleculeContainer, Domain* domain, CellProcessor* cellProcessor) :
		_decompositionTimer(decompositionTimer), _computationTimer(
				computationTimer), _domainDecomposition(domainDecomposition), _moleculeContainer(
				moleculeContainer), _domain(domain), _cellProcessor(cellProcessor) {
}

NonBlockingMPIHandlerBase::~NonBlockingMPIHandlerBase() {
}

void NonBlockingMPIHandlerBase::initBalanceAndExchange(bool forceRebalancing) {
	_decompositionTimer->start();

	_domainDecomposition->balanceAndExchange(forceRebalancing,
			_moleculeContainer, _domain);

	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateMoleculeCaches();

	_decompositionTimer->stop();
}

void NonBlockingMPIHandlerBase::performComputation(){
	_computationTimer->start();
	// Force calculation and other pair interaction related computations
	global_log->debug() << "Traversing pairs" << std::endl;
	_moleculeContainer->traverseCells(*_cellProcessor);
	_computationTimer->stop();
}
