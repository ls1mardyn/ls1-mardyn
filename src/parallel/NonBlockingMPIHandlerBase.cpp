/*
 * NonBlockingMPIBase.cpp
 *
 *  Created on: Apr 27, 2016
 *      Author: seckler
 */

#include "NonBlockingMPIHandlerBase.h"
#include "Domain.h"
#include "parallel/DomainDecompMPIBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "plugins/PluginBase.h"
#include "utils/Logger.h"


NonBlockingMPIHandlerBase::NonBlockingMPIHandlerBase(DomainDecompMPIBase* domainDecomposition,
													 ParticleContainer* moleculeContainer, Domain* domain,
													 CellProcessor* cellProcessor)
	: _domainDecomposition(domainDecomposition),
	  _moleculeContainer(moleculeContainer),
	  _domain(domain),
	  _cellProcessor(cellProcessor) {}

NonBlockingMPIHandlerBase::~NonBlockingMPIHandlerBase() {}

void NonBlockingMPIHandlerBase::performOverlappingTasks(bool forceRebalancing, double etime) {
	global_simulation->timers()->start("SIMULATION_DECOMPOSITION");
	// ensure that all Particles are in the right cells and exchange Particles
	global_log->debug() << "Updating container and decomposition" << std::endl;
	// The particles have moved, so the neighbourhood relations have
	// changed and have to be adjusted
	_moleculeContainer->update();
	global_simulation->timers()->stop("SIMULATION_DECOMPOSITION");

	// check if domain decomposition allows for non-blocking balance and exchange step. if it does, perform a
	// non-blocking step
	if (_domainDecomposition->queryBalanceAndExchangeNonBlocking(forceRebalancing, _moleculeContainer, _domain,
																 etime)) {
		// calls to derived class (if it exists, otherwise calls sequential version)
		initBalanceAndExchange(forceRebalancing, etime);
		performComputation();

	} else {
		global_log->debug()
			<< "falling back to sequential version, since domainDecomposition is blocking in this time step."
			<< std::endl;
		NonBlockingMPIHandlerBase::initBalanceAndExchange(forceRebalancing, etime);
		NonBlockingMPIHandlerBase::performComputation();
	}
}

void NonBlockingMPIHandlerBase::performComputation() {
	// Force calculation and other pair interaction related computations
	global_log->debug() << "Traversing pairs" << std::endl;

	global_simulation->timers()->start("SIMULATION_COMPUTATION");
	global_simulation->timers()->start("SIMULATION_FORCE_CALCULATION");
	_moleculeContainer->traverseCells(*_cellProcessor);
}

void NonBlockingMPIHandlerBase::initBalanceAndExchange(bool forceRebalancing, double etime) {
	global_simulation->timers()->start("SIMULATION_DECOMPOSITION");

	_domainDecomposition->balanceAndExchange(etime, forceRebalancing, _moleculeContainer, _domain);

	// The cache of the molecules must be updated/build after the exchange process,
	// as the cache itself isn't transferred
	_moleculeContainer->updateMoleculeCaches();

	global_simulation->timers()->stop("SIMULATION_DECOMPOSITION");
}
