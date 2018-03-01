/*
 * NeighbourCommunicationScheme.cpp
 *
 *  Created on: Sep 29, 2016
 *      Author: seckler
 */
//class NeighbourCommunicationScheme;
//class DirectNeighbourCommunicationScheme;
//class IndirectNeighbourCommunicationScheme;

#include "NeighbourCommunicationScheme.h"
#include "DomainDecompMPIBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "Domain.h"
#include "parallel/ZonalMethods/ZonalMethod.h"
#include <mpi.h>

#define PUSH_PULL_PARTNERS 0

NeighbourCommunicationScheme::NeighbourCommunicationScheme(unsigned int commDimms, ZonalMethod* zonalMethod) :
	_coversWholeDomain{false, false, false}, _commDimms(commDimms), _zonalMethod(zonalMethod){
	_neighbours.resize(this->getCommDims());
	
	_haloExportForceImportNeighbours.resize(this->getCommDims());
	_haloImportForceExportNeighbours.resize(this->getCommDims());
	_leavingExportNeighbours.resize(this->getCommDims());
	_leavingImportNeighbours.resize(this->getCommDims());
	
}

NeighbourCommunicationScheme::~NeighbourCommunicationScheme() {
	delete _zonalMethod;
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
	
#if PUSH_PULL_PARTNERS
	if(msgType == LEAVING_AND_HALO_COPIES) {
		msgType = LEAVING_ONLY;
		initExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);
		finalizeExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);
		
		msgType = HALO_COPIES;
	}
#endif
	initExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);

	finalizeExchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, domainDecomp);

}

void DirectNeighbourCommunicationScheme::doDirectFallBackExchange(const std::vector<HaloRegion>& haloRegions,
		MessageType msgType, DomainDecompMPIBase* domainDecomp, ParticleContainer*& moleculeContainer) { // Only Export?
	
#if PUSH_PULL_PARTNERS
	global_log->info() << "select call - initExchangeMoleculesMPI" << endl;
	selectNeighbours(msgType, false /* export */);
#endif
	
	for (const HaloRegion& haloRegion : haloRegions) {
		bool isinownprocess = true;
		for (int d = 0; d < 3; d++) {
			if (haloRegion.offset[d] && !_coversWholeDomain[d]) {
				isinownprocess = false;
			}
		}
		if (!isinownprocess) {
			continue;
		}
		// use the sequential version
		switch (msgType) {
		case LEAVING_AND_HALO_COPIES:
			// this should not be called!
			assert(false);
			break;
		case LEAVING_ONLY:
			domainDecomp->DomainDecompBase::handleDomainLeavingParticlesDirect(haloRegion, moleculeContainer);
			break;
		case HALO_COPIES:
			domainDecomp->DomainDecompBase::populateHaloLayerWithCopiesDirect(haloRegion, moleculeContainer);
			break;
		case FORCES:
			domainDecomp->DomainDecompBase::handleForceExchangeDirect(haloRegion, moleculeContainer);
			break;
		}
	}
}

void DirectNeighbourCommunicationScheme::initExchangeMoleculesMPI(ParticleContainer* moleculeContainer,
		Domain* domain, MessageType msgType, bool /*removeRecvDuplicates*/, DomainDecompMPIBase* domainDecomp) {
	// first use sequential version, if _coversWholeDomain
//	for (unsigned int d = 0; d < 3; d++) {
//		if (_coversWholeDomain[d]) {
//			// use the sequential version
//			switch (msgType) {
//			case LEAVING_AND_HALO_COPIES:
//				domainDecomp->DomainDecompBase::handleDomainLeavingParticles(d, moleculeContainer);
//				domainDecomp->DomainDecompBase::populateHaloLayerWithCopies(d, moleculeContainer);
//				break;
//			case LEAVING_ONLY:
//				domainDecomp->DomainDecompBase::handleDomainLeavingParticles(d, moleculeContainer);
//				break;
//			case HALO_COPIES:
//				domainDecomp->DomainDecompBase::populateHaloLayerWithCopies(d, moleculeContainer);
//				break;
//			case FORCES:
//				domainDecomp->DomainDecompBase::handleForceExchange(d, moleculeContainer);
//				break;
//			}
//		}
//	}
	// We mimic the direct neighbour communication also for the sequential case, otherwise things are copied multiple times or might be forgotten.
	// We have to check each direction.
	
#if PUSH_PULL_PARTNERS
	global_log->info() << "select call - initExchangeMoleculesMPI" << endl;
	selectNeighbours(msgType, false /* export */);
#endif
	
	
	double rmin[DIMgeom]; // lower corner
	double rmax[DIMgeom]; // higher corner

	for (int d = 0; d < DIMgeom; d++) {
		rmin[d] = domainDecomp->getBoundingBoxMin(d, domain);
		rmax[d] = domainDecomp->getBoundingBoxMax(d, domain);
	}
	HaloRegion ownRegion = { rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0 , global_simulation->getcutoffRadius()};
	std::vector<HaloRegion> haloRegions;
	switch (msgType) {
	case LEAVING_AND_HALO_COPIES:
		haloRegions = _zonalMethod->getLeavingExportRegions(ownRegion, global_simulation->getcutoffRadius(), _coversWholeDomain);
		doDirectFallBackExchange(haloRegions, LEAVING_ONLY, domainDecomp, moleculeContainer);
		haloRegions = _zonalMethod->getHaloExportForceImportRegions(ownRegion, global_simulation->getcutoffRadius(), _coversWholeDomain);
		doDirectFallBackExchange(haloRegions, HALO_COPIES, domainDecomp, moleculeContainer);
		break;
	case LEAVING_ONLY:
		haloRegions = _zonalMethod->getLeavingExportRegions(ownRegion, global_simulation->getcutoffRadius(), _coversWholeDomain);
		doDirectFallBackExchange(haloRegions, msgType, domainDecomp, moleculeContainer);
		break;
	case HALO_COPIES:
		haloRegions = _zonalMethod->getHaloExportForceImportRegions(ownRegion, global_simulation->getcutoffRadius(), _coversWholeDomain);
		doDirectFallBackExchange(haloRegions, msgType, domainDecomp, moleculeContainer);
		break;
	case FORCES:
		haloRegions = _zonalMethod->getHaloImportForceExportRegions(ownRegion, global_simulation->getcutoffRadius(), _coversWholeDomain);
		doDirectFallBackExchange(haloRegions, msgType, domainDecomp, moleculeContainer);
		break;
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
		Domain* /*domain*/, MessageType msgType, bool removeRecvDuplicates, DomainDecompMPIBase* domainDecomp) {

	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	cout << "enter finalizeExchangeMoleculesMPI on: " << my_rank << endl;
	
	// msg type is fixed by the fuction call, but this needs to be done for both import and export
	int numNeighbours;
	
#if PUSH_PULL_PARTNERS
	selectNeighbours(msgType, false /* export */);
	const int numExportNeighbours = _neighbours[0].size();
	selectNeighbours(msgType, true /* import */); // current _neighbours is import
	const int numImportNeighbours = _neighbours[0].size(); 
#else
	numNeighbours = _neighbours[0].size();
#endif
	
	
	// the following implements a non-blocking recv scheme, which overlaps unpacking of
	// messages with waiting for other messages to arrive
	bool allDone = false;
	double startTime = MPI_Wtime();
	
#if PUSH_PULL_PARTNERS
	numNeighbours = numImportNeighbours;
#endif

	// for 1-stage: if there is at least one neighbour with the same rank as the sending rank, make sure to remove received duplicates!
	for (int i = 0; i < numNeighbours; i++) {
		removeRecvDuplicates |= (domainDecomp->getRank() == _neighbours[0][i].getRank());
	}
	
#if PUSH_PULL_PARTNERS
	selectNeighbours(msgType, false /* export */); // last selected is export
	numNeighbours = numExportNeighbours;
	
	for (int i = 0; i < numNeighbours; i++) {
		removeRecvDuplicates |= (domainDecomp->getRank() == _neighbours[0][i].getRank());
	}
#endif
	
	
	cout << "finished removeRecvDuplicates loop on: " << my_rank << endl;

	double waitCounter = 50.0;
	double deadlockTimeOut = 360.0;
	global_log->set_mpi_output_all();
	while (not allDone) {
		allDone = true;
		
		// cout << "while loop top on: " << my_rank << endl;

		// "kickstart" processing of all Isend requests
		for (int i = 0; i < numNeighbours; ++i) { // export required (still selected)
			if (domainDecomp->getRank() != _neighbours[0][i].getRank()){
				allDone &= _neighbours[0][i].testSend(); // THIS CAUSES A SEG-FAULT
			}
		}
		
		//cout << "for 1 on: " << my_rank << endl;
		
#if PUSH_PULL_PARTNERS
		selectNeighbours(msgType, true /* import */);
		numNeighbours = numImportNeighbours;
#endif

		// get the counts and issue the Irecv-s
		for (int i = 0; i < numNeighbours; ++i) { // import required
			if (domainDecomp->getRank() != _neighbours[0][i].getRank()){
				allDone &= _neighbours[0][i].iprobeCount(domainDecomp->getCommunicator(),
					domainDecomp->getMPIParticleType());
			}

		}
		
		//cout << "for 2 on: " << my_rank << endl;

		// unpack molecules
		for (int i = 0; i < numNeighbours; ++i) { // import required (still selected)
			if (domainDecomp->getRank() != _neighbours[0][i].getRank()){
					allDone &= _neighbours[0][i].testRecv(moleculeContainer, removeRecvDuplicates, msgType==FORCES);
			}
		}
		
		//cout << "for 3 on: " << my_rank << endl;

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			global_log->warning()
					<< "DirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock warning: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << waitCounter << " seconds"
					<< std::endl;
			waitCounter += 5.0;
			for (int i = 0; i < numNeighbours; ++i) { // import required (still selected)
				if (domainDecomp->getRank() != _neighbours[0][i].getRank())
					_neighbours[0][i].deadlockDiagnosticSendRecv();
			}
			
#if PUSH_PULL_PARTNERS
			selectNeighbours(msgType, false /* export */);
			numNeighbours = numExportNeighbours;
			
			for (int i = 0; i < numNeighbours; ++i) { // import required (still selected)
				if (domainDecomp->getRank() != _neighbours[0][i].getRank())
					_neighbours[0][i].deadlockDiagnosticSendRecv();
			}
#endif
		}
		
		//cout << "for 4 on: " << my_rank << endl;

		if (waitingTime > deadlockTimeOut) {
			global_log->error()
					<< "DirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock error: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << deadlockTimeOut << " seconds"
					<< std::endl;
			for (int i = 0; i < numNeighbours; ++i) { // export required (still selected)
				if (domainDecomp->getRank() != _neighbours[0][i].getRank())
					_neighbours[0][i].deadlockDiagnosticSendRecv();
			}
			
#if PUSH_PULL_PARTNERS
			selectNeighbours(msgType, true /* import */);
			numNeighbours = numImportNeighbours;
			
			for (int i = 0; i < numNeighbours; ++i) { // export required (still selected)
				if (domainDecomp->getRank() != _neighbours[0][i].getRank())
					_neighbours[0][i].deadlockDiagnosticSendRecv();
			}
#endif
			
			Simulation::exit(457);
		}
		
		//cout << "if letztes on: " << my_rank << endl;

	} // while not allDone
	
	cout << "exit finalizeExchangeMoleculesMPI on: " << my_rank << endl;
	
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

#if PUSH_PULL_PARTNERS

void NeighbourCommunicationScheme::selectNeighbours(MessageType msgType, bool import) {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	switch(msgType) {
		case LEAVING_ONLY:
			// leavingImport / leavingExport
			global_log->info() << "selecting Neighbours LEAVING_ONLY" << endl;
			if(import) _neighbours = _leavingImportNeighbours;
			else _neighbours = _leavingExportNeighbours;
			break;
		case HALO_COPIES: 
			// haloImport / haloExport
			global_log->info() << "selecting Neighbours HALO_COPIES" << endl;
			if(import) _neighbours = _haloImportForceExportNeighbours;
			else _neighbours = _haloExportForceImportNeighbours;
			break;
		case FORCES: 
			// forceImport / forceExport
			global_log->info() << "selecting Neighbours FORCES" << endl;
			if(import) _neighbours = _haloExportForceImportNeighbours;
			else _neighbours = _haloImportForceExportNeighbours;
			break;
	}
	
	cout << "exit on: " << my_rank << endl;
}

void DirectNeighbourCommunicationScheme::shiftIfNeccessary(double *domainLength, HaloRegion *region, double *shiftArray) { // IS THIS CORRECT?
	for(int i = 0; i < 3; i++) // calculating shift 
		if(region->rmin[i] > domainLength[i])
			shiftArray[i] = -domainLength[i];
	
	for(int i = 0; i < 3; i++) // calculating shift
		if(region->rmax[i] <= 0)
			shiftArray[i] = domainLength[i];

	for(int i = 0; i < 3; i++) { // applying shift
		region->rmax[i] += shiftArray[i];
		region->rmin[i] += shiftArray[i];
	}
}

void DirectNeighbourCommunicationScheme::overlap(HaloRegion *myRegion, HaloRegion *inQuestion) { 
	/*
	 * m = myRegion, q = inQuestion, o = overlap
	 * i)  m.max < q.max ?
	 * ii) m.min < q.min ?
	 *
	 * i) | ii) | Operation
	 * -------------------------------------------
	 *  0 |  0  | o.max = q.max and o.min = m.min 
	 *  0 |  1  | o.max = q.max and o.min = q.min
	 *  1 |  0  | o.max = m.max and o.min = m.min
	 *  1 |  1  | o.max = m.max and o.min = q.min
	 * 
	 */
	HaloRegion overlap;
	
	
	for(int i = 0; i < 3; i++) {
		if(myRegion->rmax[i] < inQuestion->rmax[i]) { // 1
			if(myRegion->rmin[i] < inQuestion->rmin[i]) { // 1 1
				memcpy(overlap.rmax, myRegion->rmax, sizeof(double) * 3);
				memcpy(overlap.rmin, inQuestion->rmin, sizeof(double) * 3);
			} else { // 1 0
				memcpy(overlap.rmax, myRegion->rmax, sizeof(double) * 3);
				memcpy(overlap.rmin, myRegion->rmin, sizeof(double) * 3);
			}
		} else { // 0
			if(myRegion->rmin[i] < inQuestion->rmin[i]) { // 0 1
				memcpy(overlap.rmax, inQuestion->rmax, sizeof(double) * 3);
				memcpy(overlap.rmin, inQuestion->rmin, sizeof(double) * 3);
			} else { // 0 0
				memcpy(overlap.rmax, inQuestion->rmax, sizeof(double) * 3);
				memcpy(overlap.rmin, myRegion->rmin, sizeof(double) * 3);
			}
		}
	}
	
	// adjust width and offset?
	memcpy(inQuestion->rmax, overlap.rmax, sizeof(double) * 3);
	memcpy(inQuestion->rmin, overlap.rmin, sizeof(double) * 3);
}

bool DirectNeighbourCommunicationScheme::iOwnThis(HaloRegion* myRegion, HaloRegion* inQuestion) { // IS THIS CORRECT?
	return myRegion->rmax[0] > inQuestion->rmin[0] && myRegion->rmax[1] > inQuestion->rmin[1]
			&& myRegion->rmax[2] > inQuestion->rmin[2] && myRegion->rmin[0] <= inQuestion->rmax[0]
			&& myRegion->rmin[1] <= inQuestion->rmax[1] && myRegion->rmin[2] <= inQuestion->rmax[2];
	// myRegion->rmax > inQuestion->rmin
	// && myRegion->rmin <= inQuestion->rmax
}

/*
 * Start listening for regions of every process.
 * Send the process specific regions to every process. Tuple of process ID and region.
 * 
 * Wait till the listening got all the regions.
 * 
 * Determine who wants regions you provide.
 * --> The whos are your neighbours 
 * 
 * How many neighbours do I need to listen for? 
 * Start Listening for that amount of neighbours.
 * 
 * Send your Information to the processes who want you.
 * 
 */
void DirectNeighbourCommunicationScheme::aquireNeighbours(Domain *domain, HaloRegion *myRegion, std::vector<HaloRegion>& desiredRegions, std::vector<CommunicationPartner>& partners01, std::vector<CommunicationPartner>& partners02) {
	int my_rank; // my rank
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int num_incoming; // the number of processes in MPI_COMM_WORLD
	MPI_Comm_size(MPI_COMM_WORLD, &num_incoming);
	
	int num_regions = desiredRegions.size(); // the number of regions I would like to aquire from other processes
	
	cout << "numregions: " << num_regions << " on: " << my_rank << endl;
	
	// tell the other processes how much you are going to send
	int num_bytes_send =  sizeof(int) * 2 + (sizeof(double) * 3 + sizeof(double) * 3 + sizeof(int) * 3 + sizeof(double) * 1) * num_regions; // how many bytes am I going to send to all the other processes?
	std::vector<int> num_bytes_receive_vec(num_incoming, 0); // vector of number of bytes I am going to receive
	//MPI_Allreduce(&num_bytes_send, &num_bytes_receive, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allgather(&num_bytes_send, 1, MPI_INT, num_bytes_receive_vec.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	cout << "num_bytes_send: " << num_bytes_send << " on: " << my_rank << endl;
	cout << "exchanged number of bytes in big region exchange on rank: " << my_rank << endl;
	
	// create byte buffer
	std::vector<unsigned char> outgoing(num_bytes_send); // outgoing byte buffer
	int i = 0;
	int p = 0;
	
	// msg format: rank | number_of_regions | region_01 | region_02 | ...
	
	memcpy(outgoing.data() + i, &my_rank, sizeof(int));
	i += sizeof(int);
	memcpy(outgoing.data() + i, &num_regions, sizeof(int));
	i += sizeof(int);
	
	for(auto &region : desiredRegions) { // filling up the outgoing byte buffer
		memcpy(outgoing.data() + i, region.rmin, sizeof(double) * 3);
		i += sizeof(double) * 3;
		memcpy(outgoing.data() + i, region.rmax, sizeof(double) * 3);
		i += sizeof(double) * 3;
		memcpy(outgoing.data() + i, region.offset, sizeof(int) * 3);
		i += sizeof(int) * 3;
		memcpy(outgoing.data() + i, &region.width, sizeof(double));
		i += sizeof(double);
	}
	
	int num_bytes_receive = 0;
	std::vector<int> num_bytes_displacements(num_incoming, 0); // vector of number of bytes I am going to receive
	for (int j = 0; j < num_incoming; j++) {
		num_bytes_displacements[j] = num_bytes_receive;
		num_bytes_receive += num_bytes_receive_vec[j];
	}
	
	std::vector<unsigned char> incoming(num_bytes_receive); // the incoming byte buffer
	
	// send your regions
	//MPI_Allgather(&outgoing, num_bytes_send, MPI_BYTE, &incoming, num_bytes_receive, MPI_BYTE, MPI_COMM_WORLD);
	MPI_Allgatherv(outgoing.data(), num_bytes_send, MPI_BYTE, incoming.data(), num_bytes_receive_vec.data(), num_bytes_displacements.data(), MPI_BYTE, MPI_COMM_WORLD);
	
	cout << "exchanged desired regions on rank: " << my_rank << endl;
	
	std::vector<int> candidates(num_incoming, 0); // outgoing row
	std::vector<int> rec_information(num_incoming, 0); // how many bytes does each process expect?
	int bytesOneRegion = sizeof(double) * 3 + sizeof(double) * 3 + sizeof(int) * 3 + sizeof(double) + sizeof(double) * 3;
	std::vector<std::vector <unsigned char*>> sendingList (num_incoming); // the regions I own and want to send
	std::vector<CommunicationPartner> comm_partners02;
	
	i = 0;
	while(i != num_bytes_receive) {
		cout << i << " " << num_bytes_receive << " " << sizeof(double) << " " << sizeof(int) << "on: " << my_rank << endl;
		
		int rank;
		int regions;
		
		
		
		memcpy(&rank, incoming.data() + i, sizeof(int));
		i += sizeof(int); // 4
		memcpy(&regions, incoming.data() + i, sizeof(int));
		i += sizeof(int); // 4
		
		cout << regions << " from: " << rank << " on: " << my_rank << endl; // 26 on 1 and 0 from 1 and 0 
		
		
		for(int j = 0; j < regions; j++) {
			HaloRegion region;
			memcpy(region.rmin, incoming.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3; // 24
			memcpy(region.rmax, incoming.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3; // 24
			memcpy(region.offset, incoming.data() + i, sizeof(int) * 3);
			i += sizeof(int) * 3; // 12
			memcpy(&region.width, incoming.data() + i, sizeof(double));
			i += sizeof(double); // 4
			
			
			// msg format one region: rmin | rmax | offset | width | shift
			std::vector<double> shift(3, 0); 
			double domainLength[3] = { domain->getGlobalLength(0), domain->getGlobalLength(1), domain->getGlobalLength(2) }; // better for testing
			shiftIfNeccessary(domainLength, &region, shift.data()); 
			
			if(rank != my_rank && iOwnThis(myRegion, &region)) { 
				candidates[rank]++; // this is a region I will send to rank
				
				// cout << "candidates[rank]: " << candidates[rank] << " on: " << my_rank << endl; // 26 on both
				
				overlap(myRegion, &region); // different shift for the overlap?
				
				// make a note in partners02 - don't forget to squeeze partners02
				bool enlarged[3][2] = {{ false }};
				CommunicationPartner myNewNeighbour(rank, region.rmin, region.rmax, region.rmin, region.rmax, shift.data(), region.offset, enlarged);
				comm_partners02.push_back(myNewNeighbour);

				
				for(int k = 0; k < 3; k++) { // shift back
					region.rmax[k] -= shift[k];
					region.rmin[k] -= shift[k];
				}
				
				unsigned char* singleRegion = new unsigned char[bytesOneRegion];
							
				p = 0;
				memcpy(singleRegion + p, region.rmin, sizeof(double) * 3);
				p += sizeof(double) * 3;
				memcpy(singleRegion + p, region.rmax, sizeof(double) * 3);
				p += sizeof(double) * 3;
				memcpy(singleRegion + p, region.offset, sizeof(int) * 3);
				p += sizeof(int) * 3;
				memcpy(singleRegion + p, &region.width, sizeof(double));
				p += sizeof(double);
				memcpy(singleRegion + p, shift.data(), sizeof(double) * 3);
				p += sizeof(double) * 3;
				
				
				sendingList[rank].push_back(singleRegion); 
			}
		}
	}
	
	
	// squeeze here
	if(comm_partners02.size() > 0) {
		std::vector<CommunicationPartner> squeezed = squeezePartners(comm_partners02);
		partners02.insert(partners02.end(), squeezed.begin(), squeezed.end());
	}
	// assign here
	cout << "HERE" << endl;
	
	std::vector<unsigned char *> merged (num_incoming); // Merge each list of char arrays into one char array
	for(int j = 0; j < num_incoming; j++) {
		if(candidates[j]  > 0) {
			unsigned char* mergedRegions = new unsigned char[candidates[j] * bytesOneRegion];

			for(int k = 0; k < candidates[j]; k++) {
				memcpy(mergedRegions + k * bytesOneRegion, sendingList[j][k], bytesOneRegion);
			}
			
			merged[j] = mergedRegions;
		}
	}
	
	// delete sendingList

	for(auto one : sendingList){
		for (auto two : one){
			delete[] two;
		}
	}

	
	// The problem is, that I do not know how many regions I am going to receive from each process
	/* e.g.: 4x4
	 *
	 * row := sender
	 * column := receiver
	 *   
	 *    0 1 2 3 4
	 *   -----------
	 * 0 | | |3| | | 
	 *   -----------
	 * 1 | | | |2| |
	 *   -----------
	 * 2 | | | | | |
	 *   -----------
	 * 3 | | | | | |
	 *   -----------
	 * 
	 * Each process has a horizontal vector, where it marks how many regions it is going to send to another process.
	 * As allToAll is too slow or apparently too unpredictable regarding runtime,
	 * AllReduce is used in combination with probe.
	 * 
	 */
	
	cout << "candidates: ";
	for(int j = 0; j < num_incoming; j++) {
		cout << candidates[j] << " ";
	}
	cout << endl;
	MPI_Allreduce(candidates.data(), rec_information.data(), num_incoming,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	cout << "rec_information[my_rank]: " << rec_information[my_rank] << " on: " << my_rank << endl;
	
	// all the information for the final information exchange has been collected -> final exchange
	
	std::vector<MPI_Request> requests(num_incoming, MPI_REQUEST_NULL);
	MPI_Status probe_status;
	MPI_Status rec_status;
	
	// sending (non blocking)
	for(int j = 0; j < num_incoming; j++) {
		if(candidates[j] > 0) {
			MPI_Isend(merged[j], candidates[j] * bytesOneRegion, MPI_BYTE, j, 1, MPI_COMM_WORLD, &requests[j]); // tag is one
		}
	}
	
	//MPI_Barrier(MPI_COMM_WORLD);
	
	cout << "sent the neighbours on rank: " << my_rank << endl;
	
	std::vector<CommunicationPartner> comm_partners01; // the communication partners
	
	// receive data (blocking)
	int byte_counter = 0;
	
	while(byte_counter < rec_information[my_rank] * bytesOneRegion) { // TODO: I do not know yet, how many regions I received

		// MPI_PROBE
		MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &probe_status);
		// interpret probe
		int source = probe_status.MPI_SOURCE;
		int bytes;
		MPI_Get_count(&probe_status, MPI_BYTE, &bytes);
		
		cout << "received neighbour bytes: " << bytes << " on: " << my_rank << endl; // only one process seems to receive 2392 bytes
		
		byte_counter += bytes;
		// create buffer
		std::vector<unsigned char> raw_neighbours(bytes);
		MPI_Recv(raw_neighbours.data(), bytes, MPI_BYTE, source, 1, MPI_COMM_WORLD, &rec_status);
		// Interpret Buffer and add neighbours
		for(int k = 0; k < (bytes / bytesOneRegion); k++) { // number of regions from this process
			HaloRegion region;
			double shift[3];
			i = k * bytesOneRegion;

			memcpy(region.rmin, raw_neighbours.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3;
			memcpy(region.rmax, raw_neighbours.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3;
			memcpy(region.offset, raw_neighbours.data() + i, sizeof(int) * 3);
			i += sizeof(int) * 3;
			memcpy(&region.width, raw_neighbours.data() + i, sizeof(double));
			i += sizeof(double);

			memcpy(shift, raw_neighbours.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3;

			bool enlarged[3][2] = {{ false }};
			
			// CommunicationPartner(const int r, const double hLo[3], const double hHi[3], const double bLo[3], const double bHi[3], const double sh[3], const int offset[3], const bool enlarged[3][2]) {
			CommunicationPartner myNewNeighbour(source, region.rmin, region.rmax, region.rmin, region.rmax, shift, region.offset, enlarged); // DO NOT KNOW ABOUT THE 0s
			comm_partners01.push_back(myNewNeighbour);
			
		}
	}
	


	cout << "received the neighbours on rank: " << my_rank << endl;
	
	for(int j = 0; j < num_incoming; j++) {
		if(candidates[j] > 0) 
			MPI_Wait(&requests[j], MPI_STATUS_IGNORE);
	}

	cout << "freed the requests on rank: " << my_rank << endl;
	for (auto two : merged) {
		delete[] two;
	}
	
	cout << "number neighbousqueezePartrs: " << comm_partners01.size() << " on: " << my_rank << endl;
	if(comm_partners01.size() > 0) {
		std:vector<CommunicationPartner> squeezed = squeezePartners(comm_partners01);
		partners01.insert(partners01.end(), squeezed.begin(), squeezed.end());
		cout << "FINAL NUMBER OF NEIGHBOURS: " << partners01.size() << " on: " << my_rank << endl;
	}
	cout << "exit aquire on rank: " << my_rank << endl;


	MPI_Barrier(MPI_COMM_WORLD);
}

#endif

void DirectNeighbourCommunicationScheme::initCommunicationPartners(double cutoffRadius, Domain * domain,
		DomainDecompMPIBase* domainDecomp) {

// corners of the process-specific domain
	cout << "reached initCommunicationPartners" << endl;
	
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
	
	
	HaloRegion ownRegion = { rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0 , cutoffRadius};

#if PUSH_PULL_PARTNERS
	// halo/force regions
	std::vector<HaloRegion> haloOrForceRegions = _zonalMethod->getHaloImportForceExportRegions(ownRegion, cutoffRadius, _coversWholeDomain);
	std::vector<HaloRegion> leavingRegions = _zonalMethod->getLeavingExportRegions(ownRegion, cutoffRadius, _coversWholeDomain);
	
	// assuming p1 sends regions to p2
	aquireNeighbours(domain, &ownRegion, haloOrForceRegions, _haloImportForceExportNeighbours[0], _haloExportForceImportNeighbours[0]); // p1 notes reply, p2 notes owned as haloExportForceImport
	aquireNeighbours(domain, &ownRegion, leavingRegions, _leavingExportNeighbours[0], _leavingImportNeighbours[0]); // p1 notes reply, p2 notes owned as leaving import
	
#else
	
	std::vector<HaloRegion> haloRegions = _zonalMethod->getLeavingExportRegions(ownRegion, cutoffRadius, _coversWholeDomain);
	std::vector<CommunicationPartner> commPartners;
	for (HaloRegion haloRegion : haloRegions) {
		auto newCommPartners = domainDecomp->getNeighboursFromHaloRegion(domain, haloRegion, cutoffRadius);
		commPartners.insert(commPartners.end(), newCommPartners.begin(), newCommPartners.end());
	}
	_fullShellNeighbours = commPartners;
	//we could squeeze the fullShellNeighbours if we would want to (might however screw up FMM)
	_neighbours[0] = squeezePartners(commPartners);
#endif

}

void IndirectNeighbourCommunicationScheme::initExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer,
		Domain* /*domain*/, MessageType msgType, bool /*removeRecvDuplicates*/, unsigned short d,
		DomainDecompMPIBase* domainDecomp) {
	
	global_log->info() << "select call - initExchangeMoleculesMPI1D" << endl;
	
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

		const int numNeighbours = _neighbours[d].size();

		for (int i = 0; i < numNeighbours; ++i) {
			global_log->debug() << "Rank " << domainDecomp->getRank() << " is initiating communication to" << std::endl;
			_neighbours[d][i].initSend(moleculeContainer, domainDecomp->getCommunicator(),
					domainDecomp->getMPIParticleType(), msgType);
		}

	}
}

void IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1D(ParticleContainer* moleculeContainer,
		Domain* /*domain*/, MessageType msgType, bool removeRecvDuplicates, unsigned short d,
		DomainDecompMPIBase* domainDecomp) {
	
	global_log->info() << "select call - finalizeExchangeMoleculesMPI1D" << endl;
	
	if (_coversWholeDomain[d]) {
		return;
	}

	
	const int numNeighbours = _neighbours[d].size();
	// the following implements a non-blocking recv scheme, which overlaps unpacking of
	// messages with waiting for other messages to arrive
	bool allDone = false;
	double startTime = MPI_Wtime();

	double waitCounter = 50.0;
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
			allDone &= _neighbours[d][i].testRecv(moleculeContainer, removeRecvDuplicates, msgType==FORCES);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			global_log->warning()
					<< "IndirectNeighbourCommunicationScheme::finalizeExchangeMoleculesMPI1d: Deadlock warning: Rank "
					<< domainDecomp->getRank() << " is waiting for more than " << waitCounter << " seconds"
					<< std::endl;
			waitCounter += 5.0;
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
		DomainDecompMPIBase* domainDecomp) { // if this one is used, push pull should not (at least for now) be set

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
	convert1StageTo3StageNeighbours(commPartners, _neighbours, ownRegion, cutoffRadius);
	//squeeze neighbours -> only a single send, if rightneighbour == leftneighbour
	for (unsigned int d = 0; d < _commDimms; d++) {
		_neighbours[d]= squeezePartners(_neighbours[d]);
	}
}
