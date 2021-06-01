/*
 * NonBlockingMPIBase.hpp
 *
 *  Created on: Apr 27, 2016
 *      Author: seckler
 */
#pragma once

class DomainDecompMPIBase;
class ParticleContainer;
class Domain;
class CellProcessor;

/**
 * Provides interface for generic implementations of an overlapping Communication and Computation
 * for the balance and exchange part, in conjunction with the cell traversal.
 * This class provides blocking implementations. Derived classes should both implement the doInitBalanceAndExchange
 * method, as well as the performComputation method.
 */
class NonBlockingMPIHandlerBase {
public:
	/**
	 * Constructor of the class
	 * @param domainDecomposition domainDecomposition
	 * @param moleculeContainer the molecule container (linked cells)
	 * @param domain the domain
	 * @param cellProcessor the cell processor
	 */
	NonBlockingMPIHandlerBase(DomainDecompMPIBase* domainDecomposition, ParticleContainer* moleculeContainer,
							  Domain* domain, CellProcessor* cellProcessor);

	/**
	 * Virtual destructor of the class.
	 */
	virtual ~NonBlockingMPIHandlerBase();

	/**
	 * Call this method to perform the overlapping tasks.
	 * This method will both do the balance and exchange step of the domain decomposition, as well as the cell traversal
	 * of the cell processor. If a non-blocking version is forbidden, it falls back to a sequential version. DO NOT
	 * OVERRIDE THIS METHOD!
	 *
	 * @param forceRebalancing Defines, whether a rebalancing should be forced.
	 */
	virtual void performOverlappingTasks(bool forceRebalancing, double etime) final;

protected:
	/**
	 * The initialisation of the balance and exchange step is handled here.
	 * This should normally include the sending of the molecules.
	 * Rebalancing steps can probably not be performed here, since they are not non-blocking.
	 *
	 * @param forceRebalancing Defines, whether a rebalancing should be forced.
	 */
	virtual void initBalanceAndExchange(bool forceRebalancing, double etime);

	/**
	 * Performs the cell traversal.
	 * Derived methods should first calculate cells, that are independent of transferring molecules. O the necessary
	 * molecules arrived, further computations should be performed.
	 */
	virtual void performComputation();

	DomainDecompMPIBase* _domainDecomposition;
	ParticleContainer* _moleculeContainer;
	Domain* _domain;
	CellProcessor* _cellProcessor;
};
