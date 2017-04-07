/*
 * NeighbourCommunicationScheme.cpp
 *
 *  Created on: Sep 29, 2016
 *      Author: seckler
 */

#include "NeighbourCommunicationScheme.h"
#include "DomainDecompMPIBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "FullShell.h"
#include "Domain.h"

NeighbourCommunicationScheme::NeighbourCommunicationScheme(unsigned int commDimms) :
		_commDimms(commDimms) {
	_neighbours.resize(this->getCommDims());
	for (int d = 0; d < 3; ++d) {
		_coversWholeDomain[d] = false;
	}
	_commScheme = new FullShell();
}

NeighbourCommunicationScheme::~NeighbourCommunicationScheme() {
	delete _commScheme;
}

void DirectNeighbourCommunicationScheme::prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
		DomainDecompMPIBase* domainDecomp) {
	mardyn_assert(stageNumber < getCommDims());
	initExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);
}

void DirectNeighbourCommunicationScheme::finishNonBlockingStageImpl(ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
		DomainDecompMPIBase* domainDecomp) {
	mardyn_assert(stageNumber < getCommDims());
	finalizeExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);
}

void DirectNeighbourCommunicationScheme::exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) {
	initExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);

	finalizeExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);

}

void DirectNeighbourCommunicationScheme::initExchangeMoleculesMPI(ParticleContainer* moleculeContainer,
		Domain* /*domain*/, MessageType msgType, bool /*removeRecvDuplicates*/, DomainDecompMPIBase* domainDecomp) {
	// first use sequential version, if _coversWholeDomain
	for (unsigned int d = 0; d < 3; d++) {
		if (_coversWholeDomain[d]) {
			// use the sequential version
			switch (msgType) {
			case LEAVING_AND_HALO_COPIES:
				domainDecomp->DomainDecompBase::handleDomainLeavingParticles(d, moleculeContainer);
				domainDecomp->DomainDecompBase::populateHaloLayerWithCopies(d, moleculeContainer);
				break;
			case LEAVING_ONLY:
				domainDecomp->DomainDecompBase::handleDomainLeavingParticles(d, moleculeContainer);
				break;
			case HALO_COPIES:
				domainDecomp->DomainDecompBase::populateHaloLayerWithCopies(d, moleculeContainer);
				break;
			}
		}
	}
	// 1Stage=> only _neighbours[0] exists!
	const int numNeighbours = _neighbours[0].size();
	// send only if neighbour is actually a neighbour.
	for (int i = 0; i < numNeighbours; ++i) {
		if (_neighbours[0][i].getRank() != domainDecomp->getRank()) {
			global_log->debug() << "Rank " << domainDecomp->getRank() << "is initiating communication to";
			_neighbours[0][i].initSend(moleculeContainer, domainDecomp->getCommunicator(),
					domainDecomp->getMPIParticleType(), msgType);
		}

	}

}

void DirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI(ParticleContainer* moleculeContainer,
		Domain* /*domain*/, MessageType /*msgType*/, bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) {

	const int numNeighbours = _neighbours[0].size();
	// the following implements a non-blocking recv scheme, which overlaps unpacking of
	// messages with waiting for other messages to arrive
	bool allDone = false;
	double startTime = MPI_Wtime();

	// for 1-stage: if there is at least one neighbour with the same rank as the sending rank, make sure to remove received duplicates!
	for (int i = 0; i < numNeighbours; i++) {
		removeRecvDuplicates |= (domainDecomp->getRank() == _neighbours[0][i].getRank());
	}

	double waitCounter = 1.0;
	double deadlockTimeOut = 360.0;
	global_log->set_mpi_output_all();
	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (int i = 0; i < numNeighbours; ++i) {
			if (domainDecomp->getRank() != _neighbours[0][i].getRank())
				allDone &= _neighbours[0][i].testSend();
		}

		// get the counts and issue the Irecv-s
		for (int i = 0; i < numNeighbours; ++i) {
			if (domainDecomp->getRank() != _neighbours[0][i].getRank())
				allDone &= _neighbours[0][i].iprobeCount(domainDecomp->getCommunicator(),
						domainDecomp->getMPIParticleType());
		}

		// unpack molecules
		for (int i = 0; i < numNeighbours; ++i) {
			if (domainDecomp->getRank() != _neighbours[0][i].getRank())
				allDone &= _neighbours[0][i].testRecv(moleculeContainer, removeRecvDuplicates);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			global_log->warning()
					<< "DirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock warning: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << waitCounter << " seconds"
					<< std::endl;
			waitCounter += 1.0;
			for (int i = 0; i < numNeighbours; ++i) {
				if (domainDecomp->getRank() != _neighbours[0][i].getRank())
					_neighbours[0][i].deadlockDiagnosticSendRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			global_log->error()
					<< "DirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock error: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << deadlockTimeOut << " seconds"
					<< std::endl;
			for (int i = 0; i < numNeighbours; ++i) {
				if (domainDecomp->getRank() != _neighbours[0][i].getRank())
					_neighbours[0][i].deadlockDiagnosticSendRecv();
			}
			Simulation::exit(457);
		}

	} // while not allDone
	global_log->set_mpi_output_root(0);
}

std::vector<CommunicationPartner> squeezePartners(const std::vector<CommunicationPartner>& partners) {
	std::vector<CommunicationPartner> squeezedPartners;
	std::vector<bool> used(partners.size(), false); // flag table, that describes, whether a certain comm-partner has already been added
	for (unsigned int i = 0; i < partners.size(); i++) {
		if (used[i])
			continue;  // if we already added the neighbour, don't add it again!
		int rank = partners[i].getRank();
		CommunicationPartner tmp = partners[i];
		for (unsigned int j = i + 1; j < partners.size(); j++) {
			if (partners[j].getRank() != rank)
				continue;  // only add those with same rank
			tmp.add(partners[j]);
			used[j] = true;
		}
		squeezedPartners.push_back(tmp);
	}
	return squeezedPartners;
}

void DirectNeighbourCommunicationScheme::initCommunicationPartners(double cutoffRadius, Domain * domain,
		DomainDecompMPIBase* domainDecomp) {

// corners of the process-specific domain
	double rmin[DIMgeom]; // lower corner
	double rmax[DIMgeom]; // higher corner

	for (int d = 0; d < DIMgeom; d++) {
		rmin[d] = domainDecomp->getBoundingBoxMin(d, domain);
		rmax[d] = domainDecomp->getBoundingBoxMax(d, domain);

		// TODO: this should be safe, as long as molecules don't start flying around
		// at the speed of one cutoffRadius per time step
	}

	for (unsigned int d = 0; d < _commDimms; d++) {
		_neighbours[d].clear();
	}
	HaloRegion ownRegion = { rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0 };
	std::vector<HaloRegion> haloRegions = _commScheme->getHaloRegions(ownRegion, cutoffRadius, _coversWholeDomain);
	std::vector<CommunicationPartner> commPartners;
	for (HaloRegion haloRegion : haloRegions) {
		auto newCommPartners = domainDecomp->getNeighboursFromHaloRegion(domain, haloRegion, cutoffRadius);
		commPartners.insert(std::end(commPartners), std::begin(newCommPartners), std::end(newCommPartners));
	}
	_fullShellNeighbours = commPartners;
	//we could squeeze the fullShellNeighbours if we would want to (might however screw up FMM)
	_neighbours[0] = squeezePartners(commPartners);

}

void IndirectNeighbourCommunicationScheme::initExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer,
		Domain* /*domain*/, MessageType msgType, bool /*removeRecvDuplicates*/, unsigned short d,
		DomainDecompMPIBase* domainDecomp) {
	if (_coversWholeDomain[d]) {
		// use the sequential version

		switch (msgType) {
		case LEAVING_AND_HALO_COPIES:
			domainDecomp->DomainDecompBase::handleDomainLeavingParticles(d, moleculeContainer);
			domainDecomp->DomainDecompBase::populateHaloLayerWithCopies(d, moleculeContainer);
			break;
		case LEAVING_ONLY:
			domainDecomp->DomainDecompBase::handleDomainLeavingParticles(d, moleculeContainer);
			break;
		case HALO_COPIES:
			domainDecomp->DomainDecompBase::populateHaloLayerWithCopies(d, moleculeContainer);
			break;
		}

	} else {

		const int numNeighbours = _neighbours[d].size();

		for (int i = 0; i < numNeighbours; ++i) {
			global_log->debug() << "Rank " << domainDecomp->getRank() << " is initiating communication to" << std::endl;
			_neighbours[d][i].initSend(moleculeContainer, domainDecomp->getCommunicator(),
					domainDecomp->getMPIParticleType(), msgType);
		}

	}
}

void IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer,
		Domain* /*domain*/, MessageType /*msgType*/, bool removeRecvDuplicates, unsigned short d,
		DomainDecompMPIBase* domainDecomp) {
	if (_coversWholeDomain[d]) {
		return;
	}

	const int numNeighbours = _neighbours[d].size();
	// the following implements a non-blocking recv scheme, which overlaps unpacking of
	// messages with waiting for other messages to arrive
	bool allDone = false;
	double startTime = MPI_Wtime();

	double waitCounter = 1.0;
	double deadlockTimeOut = 360.0;
	global_log->set_mpi_output_all();
	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= _neighbours[d][i].testSend();
		}

		// get the counts and issue the Irecv-s
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= _neighbours[d][i].iprobeCount(domainDecomp->getCommunicator(),
					domainDecomp->getMPIParticleType());
		}

		// unpack molecules
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= _neighbours[d][i].testRecv(moleculeContainer, removeRecvDuplicates);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			global_log->warning()
					<< "IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock warning: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << waitCounter << " seconds"
					<< std::endl;
			waitCounter += 1.0;
			for (int i = 0; i < numNeighbours; ++i) {
				_neighbours[d][i].deadlockDiagnosticSendRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			global_log->error()
					<< "IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock error: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << deadlockTimeOut << " seconds"
					<< std::endl;
			for (int i = 0; i < numNeighbours; ++i) {
				_neighbours[d][i].deadlockDiagnosticSendRecv();
			}
			Simulation::exit(457);
		}

	} // while not allDone
	global_log->set_mpi_output_root(0);
}

void IndirectNeighbourCommunicationScheme::exchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp) {

	initExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, d, domainDecomp);

	finalizeExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, d, domainDecomp);

}

void IndirectNeighbourCommunicationScheme::exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) {

	for (unsigned short d = 0; d < getCommDims(); d++) {
		exchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, d, domainDecomp);
	}

}

void IndirectNeighbourCommunicationScheme::prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
		DomainDecompMPIBase* domainDecomp) {
	mardyn_assert(stageNumber < getCommDims());
	initExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, stageNumber, domainDecomp);
}

void IndirectNeighbourCommunicationScheme::finishNonBlockingStageImpl(ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
		DomainDecompMPIBase* domainDecomp) {
	mardyn_assert(stageNumber < getCommDims());
	finalizeExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, stageNumber, domainDecomp);
}

void IndirectNeighbourCommunicationScheme::convert1StageTo3StageNeighbours(
		const std::vector<CommunicationPartner>& commPartners,
		std::vector<std::vector<CommunicationPartner>>& neighbours, HaloRegion& ownRegion, double cutoffRadius) {
	//TODO: extend for anything else than full shell
	//TODO: implement conversion of 1StageTo3StageNeighbours

	for (const CommunicationPartner& commPartner : commPartners) {
		if (!commPartner.isFaceCommunicator()) {
			continue;  // if commPartner is not a face sharing communicator, we can ignore it!
		}
		unsigned int d = commPartner.getFaceCommunicationDirection();
		neighbours[d].push_back(commPartner);
		neighbours[d].back().enlargeInOtherDirections(d, cutoffRadius); // do this more wisely if multiple neighbours exist in that direction.
	}
}

void IndirectNeighbourCommunicationScheme::initCommunicationPartners(double cutoffRadius, Domain * domain,
		DomainDecompMPIBase* domainDecomp) {

// corners of the process-specific domain
	double rmin[DIMgeom]; // lower corner
	double rmax[DIMgeom]; // higher corner

	for (int d = 0; d < DIMgeom; d++) {
		rmin[d] = domainDecomp->getBoundingBoxMin(d, domain);
		rmax[d] = domainDecomp->getBoundingBoxMax(d, domain);

		// TODO: this should be safe, as long as molecules don't start flying around
		// at the speed of one cutoffRadius per time step
	}

	for (unsigned int d = 0; d < _commDimms; d++) {
		_neighbours[d].clear();
	}
	HaloRegion ownRegion = { rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0 };
	std::vector<HaloRegion> haloRegions = _commScheme->getHaloRegions(ownRegion, cutoffRadius, _coversWholeDomain);
	std::vector<CommunicationPartner> commPartners;
	for (HaloRegion haloRegion : haloRegions) {
		auto newCommPartners = domainDecomp->getNeighboursFromHaloRegion(domain, haloRegion, cutoffRadius);
		commPartners.insert(std::end(commPartners), std::begin(newCommPartners), std::end(newCommPartners));
	}

	_fullShellNeighbours = commPartners;
	convert1StageTo3StageNeighbours(commPartners, _neighbours, ownRegion, cutoffRadius);
	//squeeze neighbours -> only a single send, if rightneighbour == leftneighbour
	for (unsigned int d = 0; d < _commDimms; d++) {
		_neighbours[d]= squeezePartners(_neighbours[d]);
	}
}
