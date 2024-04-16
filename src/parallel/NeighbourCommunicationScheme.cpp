/*
 * NeighbourCommunicationScheme.cpp
 *
 *  Created on: Sep 29, 2016
 *      Author: seckler
 */
class NeighbourCommunicationScheme;
class DirectNeighbourCommunicationScheme;
class IndirectNeighbourCommunicationScheme;

#include <mpi.h>
#include "NeighbourCommunicationScheme.h"
#include "Domain.h"
#include "DomainDecompMPIBase.h"
#include "NeighborAcquirer.h"
#include "Simulation.h"
#include "ZonalMethods/ZonalMethod.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"

NeighbourCommunicationScheme::NeighbourCommunicationScheme(
		unsigned int commDimms, ZonalMethod* zonalMethod, bool pushPull) :
		_coversWholeDomain { false, false, false }, _commDimms(commDimms), _zonalMethod(
				zonalMethod), _pushPull(pushPull) {
	if (_pushPull) {

		_haloExportForceImportNeighbours = new std::vector<
				std::vector<CommunicationPartner>>();
		_haloImportForceExportNeighbours = new std::vector<
				std::vector<CommunicationPartner>>();
		_leavingExportNeighbours = new std::vector<
				std::vector<CommunicationPartner>>();
		_leavingImportNeighbours = new std::vector<
				std::vector<CommunicationPartner>>();

		_haloExportForceImportNeighbours->resize(this->getCommDims());
		_haloImportForceExportNeighbours->resize(this->getCommDims());
		_leavingExportNeighbours->resize(this->getCommDims());
		_leavingImportNeighbours->resize(this->getCommDims());

		_neighbours = nullptr;
	} else {
		_haloExportForceImportNeighbours = nullptr;
		_haloImportForceExportNeighbours = nullptr;
		_leavingExportNeighbours = nullptr;
		_leavingImportNeighbours = nullptr;

		_neighbours = new std::vector<std::vector<CommunicationPartner>>();
		_neighbours->resize(this->getCommDims());
	}
}

NeighbourCommunicationScheme::~NeighbourCommunicationScheme() {
	if (_pushPull) {
		delete _haloExportForceImportNeighbours;
		delete _haloImportForceExportNeighbours;
		delete _leavingExportNeighbours;
		delete _leavingImportNeighbours;
	} else {
		delete _neighbours;
	}
	delete _zonalMethod;
}

void printNeigbours(std::ofstream& stream,
		std::vector<std::vector<CommunicationPartner>>& partners) {
	for (size_t i = 0; i < partners.size(); i++) {
		auto& vector = partners[i];
		stream << "neighbor dimension: " << i << std::endl;
		for (size_t j = 0; j < vector.size(); j++) {
			auto& partner = vector[j];
			stream << "Partner: " << j << std::endl;
			partner.print(stream);
		}
	}
}

void NeighbourCommunicationScheme::printCommunicationPartners(
		std::string filename) const {
	std::ofstream checkpointfilestream;
	checkpointfilestream.open(filename.c_str());

	if (_pushPull) {
		checkpointfilestream << "haloExportForceImport:" << std::endl;
		printNeigbours(checkpointfilestream, *_haloExportForceImportNeighbours);
		checkpointfilestream << "haloImportForceExportNeighbours:" << std::endl;
		printNeigbours(checkpointfilestream, *_haloImportForceExportNeighbours);
		checkpointfilestream << "leavingExportNeighbours:" << std::endl;
		printNeigbours(checkpointfilestream, *_leavingExportNeighbours);
		checkpointfilestream << "leavingImportNeighbours:" << std::endl;
		printNeigbours(checkpointfilestream, *_leavingImportNeighbours);
	} else {
		checkpointfilestream << "neighbours:" << std::endl;
		printNeigbours(checkpointfilestream, *_neighbours);
	}
	checkpointfilestream.close();
}


void DirectNeighbourCommunicationScheme::prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
		DomainDecompMPIBase* domainDecomp) {
	mardyn_assert(stageNumber < getCommDims());
	initExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp, true);
}

void DirectNeighbourCommunicationScheme::finishNonBlockingStageImpl(ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
		DomainDecompMPIBase* domainDecomp) {
	mardyn_assert(stageNumber < getCommDims());
	finalizeExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);
}

void DirectNeighbourCommunicationScheme::exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp, bool doHaloPositionCheck) {
	//if (_pushPull) {
		if (msgType == LEAVING_AND_HALO_COPIES) {
			msgType = LEAVING_ONLY;
			initExchangeMoleculesMPI(moleculeContainer, domain, msgType,
					removeRecvDuplicates, domainDecomp, doHaloPositionCheck);
			finalizeExchangeMoleculesMPI(moleculeContainer, domain, msgType,
					removeRecvDuplicates, domainDecomp);
			moleculeContainer->deleteOuterParticles();
			msgType = HALO_COPIES;
		}
	//}
	initExchangeMoleculesMPI(moleculeContainer, domain, msgType,
			removeRecvDuplicates, domainDecomp, doHaloPositionCheck);

	finalizeExchangeMoleculesMPI(moleculeContainer, domain, msgType,
			removeRecvDuplicates, domainDecomp);

}

void DirectNeighbourCommunicationScheme::doDirectFallBackExchange(
	const std::vector<HaloRegion>& haloRegions, MessageType msgType, DomainDecompMPIBase* domainDecomp,
	ParticleContainer*& moleculeContainer, std::vector<Molecule>& invalidParticles, bool doHaloPositionCheck) {

	if (_pushPull){
		selectNeighbours(msgType, false /* export */);
	}

	for (const HaloRegion& haloRegion : haloRegions) {
		bool isInOwnProcess = true;
		for (int d = 0; d < 3; d++) {
			if (haloRegion.offset[d] && !_coversWholeDomain[d]) {
				isInOwnProcess = false;
			}
		}
		if (not isInOwnProcess) {
			continue;
		}
		// use the sequential version
		switch (msgType) {
			case LEAVING_AND_HALO_COPIES:
				// this should not be called!
				mardyn_assert(false);
				break;
			case LEAVING_ONLY:
				domainDecomp->DomainDecompBase::handleDomainLeavingParticlesDirect(haloRegion, moleculeContainer,
																				   invalidParticles);
				break;
			case HALO_COPIES:
				domainDecomp->DomainDecompBase::populateHaloLayerWithCopiesDirect(haloRegion, moleculeContainer,
																				  doHaloPositionCheck);
				break;
			case FORCES:
				domainDecomp->DomainDecompBase::handleForceExchangeDirect(haloRegion, moleculeContainer);
				break;
		}
	}
}

void DirectNeighbourCommunicationScheme::initExchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain,
																  MessageType msgType, bool /*removeRecvDuplicates*/,
																  DomainDecompMPIBase* domainDecomp,
																  bool doHaloPositionCheck) {
	// We mimic the direct neighbour communication also for the sequential case, otherwise things are copied multiple
	// times or might be forgotten.
	// We have to check each direction.

	if (_pushPull) {
		selectNeighbours(msgType, false /* export */);
	}

	// halo position check needs to be done, unless it is a invalidParticle returner and it had no invalid particles.
	auto &invalidParticles = moleculeContainer->getInvalidParticlesRef();

	if (_useSequentialFallback) {
		std::array<double, DIMgeom> rmin{};  // lower corner
		std::array<double, DIMgeom> rmax{};  // higher corner

		for (int d = 0; d < DIMgeom; d++) {
			rmin[d] = domainDecomp->getBoundingBoxMin(d, domain);
			rmax[d] = domainDecomp->getBoundingBoxMax(d, domain);
		}
		HaloRegion ownRegion = {rmin[0], rmin[1], rmin[2], rmax[0], rmax[1],
								rmax[2], 0,       0,       0,       global_simulation->getcutoffRadius()};
		std::vector<HaloRegion> haloRegions;
		double* cellLength = moleculeContainer->getHaloSize();
		std::vector<Molecule> dummy{};
		switch (msgType) {
			case LEAVING_AND_HALO_COPIES:
				haloRegions = _zonalMethod->getLeavingExportRegions(ownRegion, moleculeContainer->getCutoff(),
																	_coversWholeDomain);
				doDirectFallBackExchange(haloRegions, LEAVING_ONLY, domainDecomp, moleculeContainer, invalidParticles,
										 doHaloPositionCheck);
				haloRegions = _zonalMethod->getHaloExportForceImportRegions(ownRegion, moleculeContainer->getCutoff(),
																			_coversWholeDomain, cellLength);
				doDirectFallBackExchange(haloRegions, HALO_COPIES, domainDecomp, moleculeContainer, dummy,
										 doHaloPositionCheck);
				break;
			case LEAVING_ONLY:

				haloRegions = _zonalMethod->getLeavingExportRegions(ownRegion, moleculeContainer->getCutoff(),
																	_coversWholeDomain);
				doDirectFallBackExchange(haloRegions, msgType, domainDecomp, moleculeContainer, invalidParticles,
										 doHaloPositionCheck);
				break;
			case HALO_COPIES:
				haloRegions = _zonalMethod->getHaloExportForceImportRegions(ownRegion, moleculeContainer->getCutoff(),
																			_coversWholeDomain, cellLength);
				doDirectFallBackExchange(haloRegions, msgType, domainDecomp, moleculeContainer, dummy,
										 doHaloPositionCheck);
				break;
			case FORCES:
				haloRegions = _zonalMethod->getHaloImportForceExportRegions(ownRegion, moleculeContainer->getCutoff(),
																			_coversWholeDomain, cellLength);
				doDirectFallBackExchange(haloRegions, msgType, domainDecomp, moleculeContainer, dummy,
										 doHaloPositionCheck);
				break;
		}
	}


	// 1Stage=> only _neighbours[0] exists!
	// send only if neighbour is actually a neighbour.
	for (auto & neighbor: _neighbours->operator[](0)) {
		if (not _useSequentialFallback or neighbor.getRank() != domainDecomp->getRank()) {
			Log::global_log->debug() << "Rank " << domainDecomp->getRank() << " is initiating communication (type "
								<< msgType << ") to " << neighbor.getRank() << std::endl;
			neighbor.initSend(moleculeContainer, domainDecomp->getCommunicator(),
										  domainDecomp->getMPIParticleType(), msgType, invalidParticles, true,
										  doHaloPositionCheck);
		}
	}
	if(not invalidParticles.empty()){
		Log::global_log->error_always_output() << "NeighbourCommunicationScheme: Invalid particles that should have been "
											 "sent, are still existent. They would be lost. Aborting...\n"
										  << "BoxMin: "
										  << moleculeContainer->getBoundingBoxMin(0) << ", "
										  << moleculeContainer->getBoundingBoxMin(1) << ", "
										  << moleculeContainer->getBoundingBoxMin(2) << "\n"
										  << "BoxMax: "
										  << moleculeContainer->getBoundingBoxMax(0) << ", "
										  << moleculeContainer->getBoundingBoxMax(1) << ", "
										  << moleculeContainer->getBoundingBoxMax(2) << "\n"
										  << "The particles:" << std::endl;
		for (auto& invalidParticle : invalidParticles) {
			Log::global_log->error_always_output() << invalidParticle << std::endl;
		}
		Log::global_log->error_always_output() << "The leavingExportNeighbours:" << std::endl;
		for (auto& neighbour : (*_leavingExportNeighbours)[0]) {
			std::stringstream ss;
			neighbour.print(ss);
			Log::global_log->error_always_output() << ss.str() << std::endl;
		}
		Simulation::exit(544);
	}
}

void DirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI(ParticleContainer* moleculeContainer,
		Domain* /*domain*/, MessageType msgType, bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) {
	// msg type is fixed by the fuction call, but this needs to be done for both import and export
	int numNeighbours = 0;

	int numExportNeighbours = 0;
	int numImportNeighbours = 0;
	if (_pushPull) {
		selectNeighbours(msgType, false /* export */);
		numExportNeighbours = (*_neighbours)[0].size();
		selectNeighbours(msgType, true /* import */);  // current _neighbours is import
		numImportNeighbours = (*_neighbours)[0].size();
	} else {
		numNeighbours = (*_neighbours)[0].size();
	}

	// the following implements a non-blocking recv scheme, which overlaps unpacking of
	// messages with waiting for other messages to arrive
	bool allDone = false;
	double startTime = MPI_Wtime();
	if (_pushPull) {
		numNeighbours = numImportNeighbours;
	}

	// for 1-stage: if there is at least one neighbour with the same rank as the sending rank, make sure to remove
	// received duplicates!
	for (int i = 0; i < numNeighbours; i++) {
		removeRecvDuplicates |= (domainDecomp->getRank() == (*_neighbours)[0][i].getRank());
	}

	// local helper function to apply f to all real neighbours
	auto forAllRealNeighbors = [&](auto&& f) {
		for (auto& neighbor : (*_neighbours)[0]) {
			if (not _useSequentialFallback or domainDecomp->getRank() != neighbor.getRank()) {
				f(neighbor);
			}
		}
	};

	forAllRealNeighbors([](auto& neighbor) {
		// reset receive status
		neighbor.resetReceive();
	});

	if (_pushPull) {
		selectNeighbours(msgType, false /* export */);  // last selected is export
		numNeighbours = numExportNeighbours;

		for (int i = 0; i < numNeighbours; i++) {
			removeRecvDuplicates |= (domainDecomp->getRank() == (*_neighbours)[0][i].getRank());
		}
	}

	double waitCounter = 50.0;
	double deadlockTimeOut = 360.0;
	Log::global_log->set_mpi_output_all();
	while (not allDone) {
		allDone = true;
		if (_pushPull) {
			selectNeighbours(msgType, false /* export */);  // last selected is export
			numNeighbours = numExportNeighbours;
		}
		// "kickstart" processing of all Isend requests
		forAllRealNeighbors([&](auto& neighbor) {
			// export neighbors required (still selected)
			allDone &= neighbor.testSend();
		});

		if (_pushPull) {
			selectNeighbours(msgType, true /* import */);
			numNeighbours = numImportNeighbours;
		}

		// get the counts and issue the Irecv-s
		forAllRealNeighbors([&](auto& neighbor) {
			// import neighbors required
			allDone &= neighbor.iprobeCount(domainDecomp->getCommunicator(), domainDecomp->getMPIParticleType());
		});

		// unpack molecules
		forAllRealNeighbors([&](auto& neighbor) {
			// import neighbors required (still selected)
			allDone &= neighbor.testRecv(moleculeContainer, removeRecvDuplicates, msgType == FORCES);
		});

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			Log::global_log->warning()
				<< "DirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock warning: Rank "
				<< domainDecomp->getRank() << " is waiting for more than " << waitCounter << " seconds" << std::endl;
			waitCounter += 5.0;
			forAllRealNeighbors([&](auto& neighbor) {
				// import neighbors required (still selected)
				neighbor.deadlockDiagnosticSendRecv();
			});
			if (_pushPull) {
				selectNeighbours(msgType, false /* export */);
				numNeighbours = numExportNeighbours;
				forAllRealNeighbors([&](auto& neighbor) {
					// export neighbors required
					neighbor.deadlockDiagnosticSendRecv();
				});
			}
		}

		if (waitingTime > deadlockTimeOut) {
			Log::global_log->error()
				<< "DirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock error: Rank "
				<< domainDecomp->getRank() << " is waiting for more than " << deadlockTimeOut << " seconds"
				<< std::endl;
			forAllRealNeighbors([&](auto& neighbor) {
				// export required (still selected)
				neighbor.deadlockDiagnosticSendRecv();
			});

			if (_pushPull) {
				selectNeighbours(msgType, true /* import */);
				numNeighbours = numImportNeighbours;

				forAllRealNeighbors([&](auto& neighbor) {
					// export required (still selected)
					neighbor.deadlockDiagnosticSendRecv();
				});
			}

			Simulation::exit(457);
		}

	}  // while not allDone

	Log::global_log->set_mpi_output_root(0);
}

void NeighbourCommunicationScheme::selectNeighbours(MessageType msgType, bool import) {
	switch(msgType) {
		case LEAVING_ONLY:
			// leavingImport / leavingExport
			if(import) _neighbours = _leavingImportNeighbours;
			else _neighbours = _leavingExportNeighbours;
			break;
		case HALO_COPIES:
			// haloImport / haloExport
			if(import) _neighbours = _haloImportForceExportNeighbours;
			else _neighbours = _haloExportForceImportNeighbours;
			break;
		case FORCES:
			// forceImport / forceExport
			if(import) _neighbours = _haloExportForceImportNeighbours;
			else _neighbours = _haloImportForceExportNeighbours;
			break;
		case LEAVING_AND_HALO_COPIES:
			Log::global_log->error() << "WRONG type in selectNeighbours - this should not be used for push-pull-partners "
								   "selectNeighbours method"
								<< std::endl;
			Simulation::exit(1);
			break;
	}
}


void DirectNeighbourCommunicationScheme::initCommunicationPartners(double cutoffRadius, Domain * domain,
		DomainDecompMPIBase* domainDecomp, ParticleContainer* moleculeContainer) {
	// corners of the process-specific domain
	double rmin[DIMgeom]; // lower corner
	double rmax[DIMgeom]; // higher corner

	for (int d = 0; d < DIMgeom; d++) {
		rmin[d] = domainDecomp->getBoundingBoxMin(d, domain);
		rmax[d] = domainDecomp->getBoundingBoxMax(d, domain);

		// TODO: this should be safe, as long as molecules don't start flying around
		// at the speed of one cutoffRadius per time step
	}

	if (_pushPull) {
		for (unsigned int d = 0; d < _commDimms; d++) { // why free?
			(*_haloExportForceImportNeighbours)[d].clear();
			(*_haloImportForceExportNeighbours)[d].clear();
			(*_leavingExportNeighbours)[d].clear();
			(*_leavingImportNeighbours)[d].clear();
		}
	} else {
		for (unsigned int d = 0; d < _commDimms; d++) { // why free?
			(*_neighbours)[d].clear();
		}
	}

	HaloRegion ownRegion = {rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0, cutoffRadius};

	if (_pushPull) {
		double* cellLength = moleculeContainer->getHaloSize();
		// halo/force regions
		std::vector<HaloRegion> haloOrForceRegions =
			_zonalMethod->getHaloImportForceExportRegions(ownRegion, cutoffRadius, _coversWholeDomain, cellLength);
		std::vector<HaloRegion> leavingRegions =
				_zonalMethod->getLeavingExportRegions(ownRegion, cutoffRadius,
						_coversWholeDomain);

		std::array<double, 3> globalDomainLength{domain->getGlobalLength(0), domain->getGlobalLength(1),
												 domain->getGlobalLength(2)};
		// assuming p1 sends regions to p2
		std::tie((*_haloImportForceExportNeighbours)[0], (*_haloExportForceImportNeighbours)[0]) =
			NeighborAcquirer::acquireNeighbors(globalDomainLength, &ownRegion, haloOrForceRegions,
											   domainDecomp->getCommunicator(), _useSequentialFallback);
		// p1 notes reply, p2 notes owned as haloExportForceImport
		std::tie((*_leavingExportNeighbours)[0], (*_leavingImportNeighbours)[0]) = NeighborAcquirer::acquireNeighbors(
			globalDomainLength, &ownRegion, leavingRegions, domainDecomp->getCommunicator(), _useSequentialFallback);
		// p1 notes reply, p2 notes owned as leaving import

	} else {
		std::vector<HaloRegion> haloRegions =
				_zonalMethod->getLeavingExportRegions(ownRegion, cutoffRadius,
						_coversWholeDomain);
		std::vector<CommunicationPartner> commPartners;
		for (HaloRegion haloRegion : haloRegions) {
			auto newCommPartners = domainDecomp->getNeighboursFromHaloRegion(
					domain, haloRegion, cutoffRadius);
			commPartners.insert(commPartners.end(), newCommPartners.begin(),
					newCommPartners.end());
		}
		_fullShellNeighbours = commPartners;
		//we could squeeze the fullShellNeighbours if we would want to (might however screw up FMM)
		(*_neighbours)[0] = NeighborAcquirer::squeezePartners(commPartners);
	}

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
		case FORCES:
			domainDecomp->DomainDecompBase::handleForceExchange(d, moleculeContainer);
			break;
		}

	} else {

		const int numNeighbours = (*_neighbours)[d].size();
		std::vector<Molecule> dummy;
		for (int i = 0; i < numNeighbours; ++i) {
			Log::global_log->debug() << "Rank " << domainDecomp->getRank() << " is initiating communication to" << std::endl;
			(*_neighbours)[d][i].initSend(moleculeContainer, domainDecomp->getCommunicator(),
					domainDecomp->getMPIParticleType(), msgType, dummy, false, true/*do halo position change*/);
		}

	}
}

void IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer,
		Domain* /*domain*/, MessageType msgType, bool removeRecvDuplicates, unsigned short d,
		DomainDecompMPIBase* domainDecomp) {
	if (_coversWholeDomain[d]) {
		return;
	}
	const int numNeighbours = (*_neighbours)[d].size();
	// the following implements a non-blocking recv scheme, which overlaps unpacking of
	// messages with waiting for other messages to arrive
	bool allDone = false;
	double startTime = MPI_Wtime();

	double waitCounter = 50.0;
	double deadlockTimeOut = 360.0;
	Log::global_log->set_mpi_output_all();
	for (int i = 0; i < numNeighbours; ++i) { // reset receive status
		if (domainDecomp->getRank() != (*_neighbours)[d][i].getRank()) {
			(*_neighbours)[d][i].resetReceive();
		}
	}

	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= (*_neighbours)[d][i].testSend();
		}

		// get the counts and issue the Irecv-s
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= (*_neighbours)[d][i].iprobeCount(domainDecomp->getCommunicator(),
					domainDecomp->getMPIParticleType());
		}

		// unpack molecules
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= (*_neighbours)[d][i].testRecv(moleculeContainer, removeRecvDuplicates, msgType==FORCES);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			Log::global_log->warning()
					<< "IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock warning: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << waitCounter << " seconds"
					<< std::endl;
			waitCounter += 5.0;
			for (int i = 0; i < numNeighbours; ++i) {
				(*_neighbours)[d][i].deadlockDiagnosticSendRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			Log::global_log->error()
					<< "IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock error: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << deadlockTimeOut << " seconds"
					<< std::endl;
			for (int i = 0; i < numNeighbours; ++i) {
				(*_neighbours)[d][i].deadlockDiagnosticSendRecv();
			}
			Simulation::exit(457);
		}

	} // while not allDone
	Log::global_log->set_mpi_output_root(0);
}

void IndirectNeighbourCommunicationScheme::exchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp) {

	initExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, d, domainDecomp);

	finalizeExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, d, domainDecomp);

}

void IndirectNeighbourCommunicationScheme::exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp, bool /*doHaloPositionCheck*/) {
	for (unsigned int d = 0; d < getCommDims(); d++) {
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
		DomainDecompMPIBase* domainDecomp,ParticleContainer* /*moleculeContainer*/) { // if this one is used, push pull should not (at least for now) be set

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
		(*_neighbours)[d].clear();
	}
	HaloRegion ownRegion = { rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0, cutoffRadius}; // region of the box
	// ---
	std::vector<HaloRegion> haloRegions = _zonalMethod->getLeavingExportRegions(ownRegion, cutoffRadius, _coversWholeDomain); // halo regions (outside of the box)
	// ---
	std::vector<CommunicationPartner> commPartners;
	for (HaloRegion haloRegion : haloRegions) { // determine who to communicate with - who's region is in the haloRegions vector?
		auto newCommPartners = domainDecomp->getNeighboursFromHaloRegion(domain, haloRegion, cutoffRadius);
		commPartners.insert(commPartners.end(), newCommPartners.begin(), newCommPartners.end());
	}

	_fullShellNeighbours = commPartners;
	convert1StageTo3StageNeighbours(commPartners, (*_neighbours), ownRegion, cutoffRadius);
	//squeeze neighbours -> only a single send, if rightneighbour == leftneighbour
	for (unsigned int d = 0; d < _commDimms; d++) {
		(*_neighbours)[d]= NeighborAcquirer::squeezePartners((*_neighbours)[d]);
	}
}
