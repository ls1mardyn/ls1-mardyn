/*
 * NonBlockingMPIBase.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: seckler
 */
#pragma once

class Timer;
class DomainDecompBase;
class ParticleContainer;
class Domain;
class CellProcessor;

/**
 * Provides interface for generic implementations of an overlapping Communication and Computation
 * for the balance and exchange part, in conjunction with the cell traversal.
 * This class provides blocking implementations. Derived classes should both implement the initBalanceAndExchange method,
 * as well as the performComputation method.
 */
class NonBlockingMPIHandlerBase {
public:
	/**
	 * Constructor of the class
	 * @param decompositionTimer timer for the decomposition
	 * @param computationTimer timer for the computation
	 * @param domainDecomposition domainDecomposition
	 * @param moleculeContainer the molecule container (linked cells)
	 * @param domain the domain
	 * @param cellProcessor the cell processor
	 */
	NonBlockingMPIHandlerBase(Timer* decompositionTimer,
			Timer* computationTimer, DomainDecompBase* domainDecomposition,
			ParticleContainer* moleculeContainer, Domain* domain,
			CellProcessor* cellProcessor);

	/**
	 * Virtual destructor of the class.
	 */
	virtual ~NonBlockingMPIHandlerBase();

	/**
	 * The initialisation of the balance and exchange step is handled here.
	 * This should normally include the sending of the molecules.
	 * Rebalancing steps can probably not be performed here, since they are not non-blocking.
	 *
	 * @param forceRebalancing Defines, whether a rebalancing should be forced.
	 */
	virtual void initBalanceAndExchange(bool forceRebalancing);

	/**
	 * Performs the cell traversal.
	 * Derived methods should first calculate cells, that are independent of transferring molecules. After the necessary molecules arrived,
	 * further computations should be performed.
	 */
	virtual void performComputation();

private:
	Timer* _decompositionTimer, *_computationTimer;
	DomainDecompBase* _domainDecomposition;
	ParticleContainer* _moleculeContainer;
	Domain* _domain;
	CellProcessor* _cellProcessor;
};

