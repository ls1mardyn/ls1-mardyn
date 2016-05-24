/*
 * DomainDecompBaseMPI.cpp
 *
 *  Created on: Nov 15, 2015
 *      Author: tchipevn
 */

#include "DomainDecompMPIBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"

using Log::global_log;

DomainDecompMPIBase::DomainDecompMPIBase() : _comm(MPI_COMM_WORLD) {
	for (int d = 0; d < 3; ++d) {
		_coversWholeDomain[d] = false;
	}

	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &_rank) );

	MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &_numProcs) );

	ParticleData::setMPIType(_mpiParticleType);
}

DomainDecompMPIBase::~DomainDecompMPIBase() {
	MPI_Type_free(&_mpiParticleType);

	// MPI_COMM_WORLD doesn't need to be freed, so
	// if a derived class does something with the communicator
	// then the derived class should also free it
}

unsigned DomainDecompMPIBase::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	unsigned* moldistribution = new unsigned[_numProcs];
	MPI_CHECK( MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution, 1, MPI_UNSIGNED, _comm) );
	unsigned globalN = 0;
	for (int r = 0; r < _rank; r++)
		globalN += moldistribution[r];
	unsigned localNbottom = globalN;
	globalN += moldistribution[_rank];
	unsigned localNtop = globalN;
	for (int r = _rank + 1; r < _numProcs; r++)
		globalN += moldistribution[r];
	delete[] moldistribution;
	*minrnd = (float) localNbottom / globalN;
	*maxrnd = (float) localNtop / globalN;
	return globalN;
}

void DomainDecompMPIBase::assertIntIdentity(int IX) {
	if (_rank)
		MPI_CHECK( MPI_Send(&IX, 1, MPI_INT, 0, 2 * _rank + 17, _comm) );
	else {
		int recv;
		MPI_Status s;
		for (int i = 1; i < _numProcs; i++) {
			MPI_CHECK( MPI_Recv(&recv, 1, MPI_INT, i, 2 * i + 17, _comm, &s) );
			if (recv != IX) {
				global_log->error() << "IX is " << IX << " for rank 0, but " << recv << " for rank " << i << ".\n";
				MPI_Abort(MPI_COMM_WORLD, 911);
			}
		}
		global_log->info() << "IX = " << recv << " for all " << _numProcs << " ranks.\n";
	}
}

void DomainDecompMPIBase::assertDisjunctivity(TMoleculeContainer* mm) {
	using std::map;
	using std::endl;

	Molecule* m;

	if (_rank) {
		int num_molecules = mm->getNumberOfParticles();
		unsigned long *tids;
		tids = new unsigned long[num_molecules];

		int i = 0;
		for (m = mm->begin(); m != mm->end(); m = mm->next()) {
			tids[i] = m->id();
			i++;
		}
		MPI_CHECK( MPI_Send(tids, num_molecules, MPI_UNSIGNED_LONG, 0, 2674 + _rank, _comm) );
		delete[] tids;
		global_log->info() << "Data consistency checked: for results see rank 0." << endl;
	}
	else {
		map<unsigned long, int> check;

		for (m = mm->begin(); m != mm->end(); m = mm->next())
			check[m->id()] = 0;

		MPI_Status status;
		for (int i = 1; i < _numProcs; i++) {
			int num_recv = 0;
			unsigned long *recv;
			MPI_CHECK( MPI_Probe(i, 2674 + i, _comm, &status) );
			MPI_CHECK( MPI_Get_count(&status, MPI_UNSIGNED_LONG, &num_recv) );
			recv = new unsigned long[num_recv];

			MPI_CHECK( MPI_Recv(recv, num_recv, MPI_UNSIGNED_LONG, i, 2674 + i, _comm, &status) );
			for (int j = 0; j < num_recv; j++) {
				if (check.find(recv[j]) != check.end()) {
					global_log->error() << "Ranks " << check[recv[j]] << " and " << i << " both propagate ID " << recv[j] << endl;
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				else
					check[recv[j]] = i;
			}
			delete[] recv;
		}
		global_log->info() << "Data consistency checked: No duplicate IDs detected among " << check.size() << " entries." << endl;
	}
}

void DomainDecompMPIBase::balanceAndExchangeInitNonBlocking(
		bool /*forceRebalancing*/, ParticleContainer* /*moleculeContainer*/,
		Domain* /*domain*/) {
	// for now, nothing to be done here
	// later switch between different communication schemes might go in here.
}

void DomainDecompMPIBase::prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates){
	if(DomainDecompBase::getNonBlockingStageCount() == 3){
		initExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, stageNumber);
	}
}

void DomainDecompMPIBase::finishNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates){
	if(DomainDecompBase::getNonBlockingStageCount() == 3){
		finalizeExchangeMoleculesMPI1D(moleculeContainer, domain, msgType, removeRecvDuplicates, stageNumber);
	}
}

void DomainDecompMPIBase::initExchangeMoleculesMPI1D(
		ParticleContainer* moleculeContainer, Domain* /*domain*/,
		MessageType msgType, bool /*removeRecvDuplicates*/, unsigned short d) {
	if (_coversWholeDomain[d]) {
		// use the sequential version

		switch (msgType) {
		case LEAVING_AND_HALO_COPIES:
			DomainDecompBase::handleDomainLeavingParticles(d,
					moleculeContainer);
			DomainDecompBase::populateHaloLayerWithCopies(d, moleculeContainer);
			break;
		case LEAVING_ONLY:
			DomainDecompBase::handleDomainLeavingParticles(d,
					moleculeContainer);
			break;
		case HALO_COPIES:
			DomainDecompBase::populateHaloLayerWithCopies(d, moleculeContainer);
			break;
		}
		return;
	}

	const int numNeighbours = _neighbours[d].size();

	for (int i = 0; i < numNeighbours; ++i) {
		global_log->debug() << "Rank " << _rank
				<< "is initiating communication to";
		_neighbours[d][i].initSend(moleculeContainer, _comm, _mpiParticleType,
				msgType);
	}
}

void DomainDecompMPIBase::finalizeExchangeMoleculesMPI1D(
		ParticleContainer* moleculeContainer, Domain* /*domain*/,
		MessageType msgType, bool removeRecvDuplicates, unsigned short d) {
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

	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= _neighbours[d][i].testSend();
		}

		// get the counts and issue the Irecv-s
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= _neighbours[d][i].iprobeCount(_comm, _mpiParticleType);
		}

		// unpack molecules
		for (int i = 0; i < numNeighbours; ++i) {
			allDone &= _neighbours[d][i].testRecv(moleculeContainer,
					removeRecvDuplicates);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			global_log->warning() << "Deadlock warning: Rank " << _rank
					<< " is waiting for more than " << waitCounter << " seconds"
					<< std::endl;
			waitCounter += 1.0;
			for (int i = 0; i < numNeighbours; ++i) {
				_neighbours[d][i].deadlockDiagnosticSendRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			global_log->warning() << "Deadlock error: Rank " << _rank
					<< " is waiting for more than " << deadlockTimeOut
					<< " seconds" << std::endl;
			for (int i = 0; i < numNeighbours; ++i) {
				_neighbours[d][i].deadlockDiagnosticSendRecv();
			}
			MPI_Abort(_comm, 1);
			exit(1);
		}

	} // while not allDone
}


void DomainDecompMPIBase::exchangeMoleculesMPI1D(
		ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates, unsigned short d) {

	initExchangeMoleculesMPI1D(moleculeContainer, domain, msgType,
					removeRecvDuplicates, d);

	finalizeExchangeMoleculesMPI1D(moleculeContainer, domain, msgType,
						removeRecvDuplicates, d);

}

void DomainDecompMPIBase::exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain, MessageType msgType, bool removeRecvDuplicates) {

	global_log->set_mpi_output_all();

	for (unsigned short d = 0; d < DIM; d++) {
		exchangeMoleculesMPI1D(moleculeContainer, domain, msgType,
				removeRecvDuplicates, d);
	} // for d

	global_log->set_mpi_output_root(0);
}

