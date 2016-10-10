/*
 * NeighbourCommunicationScheme.h
 *
 *  Created on: Sep 29, 2016
 *      Author: seckler
 */

#pragma once

#include "parallel/CommunicationPartner.h"

#include <vector>

class DomainDecompMPIBase;
class Domain;

class NeighbourCommunicationScheme {
public:
	/**
	 * Specifies the amount of sequential communication steps needed for the communication scheme.
	 * This is also the outer size of DomainDecompMPIBase::_neighbours
	 * @return
	 */
	unsigned int getCommDims() {
		return _commDimms;
	}
	NeighbourCommunicationScheme() = delete;
	NeighbourCommunicationScheme(unsigned int commDimms);
	virtual ~NeighbourCommunicationScheme();

	virtual void prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
			DomainDecompMPIBase* domainDecomp);

	virtual void finishNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
			unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
			DomainDecompMPIBase* domainDecomp);

	virtual void exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) = 0;

	void setCoverWholeDomain(unsigned int d, bool covers) {
		_coversWholeDomain[d] = covers;
	}

	virtual void initCommunicationPartners(double cutoffRadius, Domain * domain, DomainDecompMPIBase* domainDecomp) = 0;
protected:

	//! vector of neighbours. The first dimension should be of size getCommDims().
	std::vector<std::vector<CommunicationPartner>> _neighbours;

	//! flag, which tells whether a processor covers the whole domain along a dimension
	//! if true, we will use the methods provided by the base class for handling the
	//! respective dimension, instead of packing and unpacking messages to self
	bool _coversWholeDomain[3];

	unsigned int _commDimms;

};

class NeighbourCommunicationScheme1Stage: NeighbourCommunicationScheme {
public:
	NeighbourCommunicationScheme1Stage() :
			NeighbourCommunicationScheme(1) {
	}
	virtual ~NeighbourCommunicationScheme1Stage();
};

class NeighbourCommunicationScheme3Stage: NeighbourCommunicationScheme {
public:

	NeighbourCommunicationScheme3Stage() :
			NeighbourCommunicationScheme(3) {
	}
	virtual ~NeighbourCommunicationScheme3Stage();
	void exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp);
	virtual void initCommunicationPartners(double cutoffRadius, Domain * domain, DomainDecompMPIBase* domainDecomp);
protected:
	void initExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp);

	void finalizeExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp);
	void exchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType,
			bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp);

};
