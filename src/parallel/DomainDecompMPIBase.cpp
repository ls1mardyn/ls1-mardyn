/*
 * DomainDecompBaseMPI.cpp
 *
 *  Created on: Nov 15, 2015
 *      Author: tchipevn
 */

#include "DomainDecompMPIBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "parallel/NeighbourCommunicationScheme.h"
#include "ParticleData.h"

using Log::global_log;

DomainDecompMPIBase::DomainDecompMPIBase() :
		_comm(MPI_COMM_WORLD) {

	_neighbourCommunicationScheme = new IndirectNeighbourCommunicationScheme();
	//_neighbourCommunicationScheme = new DirectNeighbourCommunicationScheme();

	MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &_rank));

	MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &_numProcs));

	ParticleData::getMPIType(_mpiParticleType);
}

DomainDecompMPIBase::~DomainDecompMPIBase() {

	delete _neighbourCommunicationScheme;
	MPI_Type_free(&_mpiParticleType);

	// MPI_COMM_WORLD doesn't need to be freed, so
	// if a derived class does something with the communicator
	// then the derived class should also free it
}

void DomainDecompMPIBase::readXML(XMLfileUnits& xmlconfig) {
	std::string communicationScheme = "indirect";
	xmlconfig.getNodeValue("CommunicationScheme", communicationScheme);
	setCommunicationScheme(communicationScheme);
}

int DomainDecompMPIBase::getNonBlockingStageCount(){
	return _neighbourCommunicationScheme->getCommDims();
}

void DomainDecompMPIBase::setCommunicationScheme(std::string scheme){
	if(_neighbourCommunicationScheme!=nullptr){
		delete _neighbourCommunicationScheme;
	}
	if (scheme=="direct"){
		global_log->info() << "DomainDecompMPIBase: Using DirectCommunicationScheme" << std::endl;
		_neighbourCommunicationScheme = new DirectNeighbourCommunicationScheme();
	} else if(scheme=="indirect"){
		global_log->info() << "DomainDecompMPIBase: Using IndirectCommunicationScheme" << std::endl;
		_neighbourCommunicationScheme = new IndirectNeighbourCommunicationScheme();
	} else{
		global_log->error() << "DomainDecompMPIBase: invalid CommunicationScheme specified. Valid values are 'direct' and 'indirect'"
				<< std::endl;
		Simulation::exit(1);
	}
}

unsigned DomainDecompMPIBase::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	unsigned* moldistribution = new unsigned[_numProcs];
	MPI_CHECK(MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution, 1, MPI_UNSIGNED, _comm));
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
		MPI_CHECK(MPI_Send(&IX, 1, MPI_INT, 0, 2 * _rank + 17, _comm));
	else {
		int recv;
		MPI_Status s;
		for (int i = 1; i < _numProcs; i++) {
			MPI_CHECK(MPI_Recv(&recv, 1, MPI_INT, i, 2 * i + 17, _comm, &s));
			if (recv != IX) {
				global_log->error() << "IX is " << IX << " for rank 0, but " << recv << " for rank " << i << ".\n";
				MPI_Abort(MPI_COMM_WORLD, 911);
			}
		}
		global_log->info() << "IX = " << recv << " for all " << _numProcs << " ranks.\n";
	}
}

void DomainDecompMPIBase::assertDisjunctivity(TMoleculeContainer* mm) const {
	using std::map;
	using std::endl;

	if (_rank) {
		unsigned long num_molecules = mm->getNumberOfParticles();
		unsigned long *tids = new unsigned long[num_molecules];

		int i = 0;
		for (ParticleIterator m = mm->iteratorBegin(); m != mm->iteratorEnd(); ++m) {
			tids[i] = m->id();
			i++;
		}
		MPI_CHECK(MPI_Send(tids, num_molecules, MPI_UNSIGNED_LONG, 0, 2674 + _rank, _comm));
		delete[] tids;
		global_log->info() << "Data consistency checked: for results see rank 0." << endl;
	} else {
		map<unsigned long, int> check;

		for (ParticleIterator m = mm->iteratorBegin(); m != mm->iteratorEnd(); ++m)
			check[m->id()] = 0;

		MPI_Status status;
		for (int i = 1; i < _numProcs; i++) {
			int num_recv = 0;
			unsigned long *recv;
			MPI_CHECK(MPI_Probe(i, 2674 + i, _comm, &status));
			MPI_CHECK(MPI_Get_count(&status, MPI_UNSIGNED_LONG, &num_recv));
			recv = new unsigned long[num_recv];

			MPI_CHECK(MPI_Recv(recv, num_recv, MPI_UNSIGNED_LONG, i, 2674 + i, _comm, &status));
			for (int j = 0; j < num_recv; j++) {
				if (check.find(recv[j]) != check.end()) {
					global_log->error() << "Ranks " << check[recv[j]] << " and " << i << " both propagate ID "
							<< recv[j] << endl;
					MPI_Abort(MPI_COMM_WORLD, 1);
				} else
					check[recv[j]] = i;
			}
			delete[] recv;
		}
		global_log->info() << "Data consistency checked: No duplicate IDs detected among " << check.size()
				<< " entries." << endl;
	}
}

void DomainDecompMPIBase::balanceAndExchangeInitNonBlocking(bool /*forceRebalancing*/,
		ParticleContainer* /*moleculeContainer*/, Domain* /*domain*/) {
	// for now, nothing to be done here
	// later switch between different communication schemes might go in here.
}

void DomainDecompMPIBase::prepareNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates) {
	mardyn_assert(stageNumber < _neighbourCommunicationScheme->getCommDims());
	_neighbourCommunicationScheme->prepareNonBlockingStageImpl(moleculeContainer, domain, stageNumber, msgType,
			removeRecvDuplicates, this);
}

void DomainDecompMPIBase::finishNonBlockingStageImpl(ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates) {
	_neighbourCommunicationScheme->finishNonBlockingStageImpl(moleculeContainer, domain, stageNumber, msgType,
			removeRecvDuplicates, this);
}

void DomainDecompMPIBase::exchangeMoleculesMPI(ParticleContainer* moleculeContainer, Domain* domain,
		MessageType msgType, bool removeRecvDuplicates) {

	global_log->set_mpi_output_all();

	_neighbourCommunicationScheme->exchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates, this);

	global_log->set_mpi_output_root(0);
}

size_t DomainDecompMPIBase::getTotalSize() {
	return DomainDecompBase::getTotalSize() + _neighbourCommunicationScheme->getDynamicSize()
			+ _collCommunication.getDynamicSize();
}

void DomainDecompMPIBase::printSubInfo(int offset){
	std::stringstream offsetstream;
	for (int i = 0; i < offset; i++) {
		offsetstream << "\t";
	}
	global_log->info() << offsetstream.str() << "own datastructures:\t" << sizeof(DomainDecompMPIBase) / 1.e6 << " MB" << std::endl;
	global_log->info() << offsetstream.str() << "neighbourCommunicationScheme:\t\t" << _neighbourCommunicationScheme->getDynamicSize() / 1.e6 << " MB" << std::endl;
	global_log->info() << offsetstream.str() << "collective Communication:\t\t" << _collCommunication.getDynamicSize() / 1.e6 << " MB" << std::endl;

}
