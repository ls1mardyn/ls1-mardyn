/*
 * NonBlockingMPIMultiStepHandler.h
 *
 *  Created on: Apr 28, 2016
 *      Author: seckler
 */

#pragma once

#include "NonBlockingMPIHandlerBase.h"

class NonBlockingMPIMultiStepHandler : public NonBlockingMPIHandlerBase {
public:
	/**
	 * Constructor of the class
	 * @param domainDecomposition domainDecomposition
	 * @param moleculeContainer the molecule container (linked cells)
	 * @param domain the domain
	 * @param cellProcessor the cell processor
	 */
	NonBlockingMPIMultiStepHandler(DomainDecompMPIBase* domainDecomposition, ParticleContainer* moleculeContainer,
								   Domain* domain, CellProcessor* cellProcessor);

	/**
	 * Virtual destructor of the class.
	 */
	virtual ~NonBlockingMPIMultiStepHandler();

protected:
	// documentation in base class
	virtual void initBalanceAndExchange(bool forceRebalancing, double etime);

	// documentation in base class
	virtual void performComputation();
};
