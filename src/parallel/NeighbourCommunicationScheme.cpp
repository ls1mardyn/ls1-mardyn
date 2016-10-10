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

void NeighbourCommunicationScheme3Stage::initExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer,
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
			global_log->debug() << "Rank " << domainDecomp->getRank() << "is initiating communication to";
			_neighbours[d][i].initSend(moleculeContainer, domainDecomp->getCommunicator(),
					domainDecomp->getMPIParticleType(), msgType);
		}
	}
}

void NeighbourCommunicationScheme3Stage::finalizeExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer,
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
	double deadlockTimeOut = 60.0;
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
			global_log->warning() << "DomainDecompMPIBase::finalizeExchangeMoleculesMPI1d: Deadlock warning: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << waitCounter << " seconds"
					<< std::endl;
			waitCounter += 1.0;
			for (int i = 0; i < numNeighbours; ++i) {
				_neighbours[d][i].deadlockDiagnosticSendRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			global_log->error() << "DomainDecompMPIBase::finalizeExchangeMoleculesMPI1d: Deadlock error: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << deadlockTimeOut << " seconds"
					<< std::endl;
			for (int i = 0; i < numNeighbours; ++i) {
				_neighbours[d][i].deadlockDiagnosticSendRecv();
			}
			global_simulation->exit(457);
		}

	} // while not allDone
	global_log->set_mpi_output_root(0);
}

void NeighbourCommunicationScheme3Stage::exchangeMoleculesMPI1D(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, unsigned short d, DomainDecompMPIBase* domainDecomp) {

	initExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, d, domainDecomp);

	finalizeExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, d, domainDecomp);

}

void NeighbourCommunicationScheme3Stage::exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) {
	for (unsigned short d = 0; d < getCommDims(); d++) {
		exchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, d, domainDecomp);
	}
}

void NeighbourCommunicationScheme3Stage::initCommunicationPartners(double cutoffRadius, Domain * domain,
		DomainDecompMPIBase* domainDecomp) {

	// corners of the process-specific domain
	double rmin[DIMgeom]; // lower corner
	double rmax[DIMgeom]; // higher corner
	double halo_width[DIMgeom]; // width of the halo strip

	for (int d = 0; d < DIMgeom; d++) {
		rmin[d] = domainDecomp->getBoundingBoxMin(d, domain);
		rmax[d] = domainDecomp->getBoundingBoxMax(d, domain);

		// TODO: this should be safe, as long as molecules don't start flying around
		// at the speed of one cutoffRadius per time step
		halo_width[d] = cutoffRadius;
	}

	for (unsigned int d = 0; d < _commDimms; d++) {
		_neighbours[d].clear();
	}
	HaloRegion ownRegion = { rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0 };
	std::vector<HaloRegion> haloRegions = _commScheme->getHaloRegions(ownRegion, cutoffRadius, _coversWholeDomain);
	std::vector<CommunicationPartner> commPartners;
	for (HaloRegion haloRegion : haloRegions) {
		commPartners.push_back(domainDecomp->getNeighboursFromHaloRegion(haloRegion));
	}

	for (unsigned short dimension = 0; dimension < DIMgeom; dimension++) {
		if (_coversWholeDomain[dimension]) {
			// nothing to do;
			continue;
		}

		// set the ranks
		int ranks[2];

		MPI_CHECK(MPI_Cart_shift(domainDecomp->getCommunicator(), dimension, 1, &ranks[LOWER], &ranks[HIGHER]));

		// When moving a particle across a periodic boundary, the molecule position has to change.
		// These offsets specify for each dimension (x, y and z) and each direction ("left"/lower
		// neighbour and "right"/higher neighbour, how the particle coordinates have to be changed.
		// e.g. for dimension x (d=0) and a process on the left boundary of the domain, particles
		// moving to the left get the length of the whole domain added to their x-value
		double offsetLower[DIMgeom];
		double offsetHigher[DIMgeom];
		offsetLower[dimension] = 0.0;
		offsetHigher[dimension] = 0.0;

		// process on the left boundary
		if (_coords[dimension] == 0)
			offsetLower[dimension] = domain->getGlobalLength(dimension);
		// process on the right boundary
		if (_coords[dimension] == _gridSize[dimension] - 1)
			offsetHigher[dimension] = -domain->getGlobalLength(dimension);

		for (int direction = LOWER; direction <= HIGHER; direction++) {
			double regToSendLow[DIMgeom];
			double regToSendHigh[DIMgeom];

			// set the regions
			for (int i = 0; i < DIMgeom; i++) {
				regToSendLow[i] = rmin[i] - halo_width[i];
				regToSendHigh[i] = rmax[i] + halo_width[i];
			}

			double haloLow[3];
			double haloHigh[3];
			double boundaryLow[3];
			double boundaryHigh[3];

			switch (direction) {
			case LOWER:
				regToSendHigh[dimension] = rmin[dimension] + halo_width[dimension];
				for (int i = 0; i < DIMgeom; ++i) {
					haloLow[i] = regToSendLow[i];
					if (i == dimension) {
						haloHigh[i] = boundaryLow[i] = rmin[i];
					} else {
						haloHigh[i] = regToSendHigh[i];
						boundaryLow[i] = regToSendLow[i];
					}
					boundaryHigh[i] = regToSendHigh[i];
				}
				break;
			case HIGHER:
				regToSendLow[dimension] = rmax[dimension] - halo_width[dimension];
				for (int i = 0; i < DIMgeom; ++i) {
					boundaryLow[i] = regToSendLow[i];
					if (i == dimension) {
						boundaryHigh[i] = haloLow[i] = rmax[i];
					} else {
						boundaryHigh[i] = regToSendHigh[i];
						haloLow[i] = regToSendLow[i];
					}
					haloHigh[i] = regToSendHigh[i];
				}
				break;
			}

			// set the shift
			double shift[3] = { 0., 0., 0. };
			if (direction == LOWER)
				shift[dimension] = offsetLower[dimension];
			if (direction == HIGHER)
				shift[dimension] = offsetHigher[dimension];
			_neighbours[dimension].push_back(
					CommunicationPartner(ranks[direction], haloLow, haloHigh, boundaryLow, boundaryHigh, shift));
		}
	}
}
