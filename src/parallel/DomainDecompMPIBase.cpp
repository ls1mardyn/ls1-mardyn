/*
 * DomainDecompBaseMPI.cpp
 *
 *  Created on: Nov 15, 2015
 *	  Author: tchipevn
 */
#include <memory>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "DomainDecompMPIBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "parallel/NeighbourCommunicationScheme.h"
#include "ParticleData.h"
#include "Domain.h"

#include "parallel/ZonalMethods/FullShell.h"
#include "parallel/ZonalMethods/EighthShell.h"
#include "parallel/ZonalMethods/HalfShell.h"
#include "parallel/ZonalMethods/Midpoint.h"
#include "parallel/ZonalMethods/NeutralTerritory.h"
#include "parallel/CollectiveCommunication.h"
#include "parallel/CollectiveCommunicationNonBlocking.h"


DomainDecompMPIBase::DomainDecompMPIBase() : DomainDecompMPIBase(MPI_COMM_WORLD) {}

DomainDecompMPIBase::DomainDecompMPIBase(MPI_Comm comm) :
		_comm(comm) {
#ifndef MARDYN_AUTOPAS
	_neighbourCommunicationScheme = std::make_unique<IndirectNeighbourCommunicationScheme>(new FullShell());
#else
	// direct push-pull
	_neighbourCommunicationScheme = std::make_unique<DirectNeighbourCommunicationScheme>(new FullShell(), true);
#endif
	//_neighbourCommunicationScheme = new DirectNeighbourCommunicationScheme(new FullShell());

	MPI_CHECK(MPI_Comm_rank(_comm, &_rank));

	MPI_CHECK(MPI_Comm_size(_comm, &_numProcs));

	ParticleData::getMPIType(_mpiParticleType);

	_collCommunication = std::unique_ptr<CollectiveCommunicationInterface>(new CollectiveCommunication());
	//_collCommunication = std::unique_ptr<CollectiveCommunicationInterface>(new CollectiveCommunicationNonBlocking());
}

DomainDecompMPIBase::~DomainDecompMPIBase() {

	MPI_Type_free(&_mpiParticleType);

	// MPI_COMM_WORLD doesn't need to be freed, so
	// if a derived class does something with the communicator
	// then the derived class should also free it
}

void DomainDecompMPIBase::readXML(XMLfileUnits& xmlconfig) {
	// store current path
	std::string oldPath(xmlconfig.getcurrentnodepath());

#ifndef MARDYN_AUTOPAS
	std::string neighbourCommunicationScheme = "indirect";
	xmlconfig.getNodeValue("CommunicationScheme", neighbourCommunicationScheme);
#else
	std::string neighbourCommunicationScheme = "direct-pp";
#endif
	if(_forceDirectPP){
		Log::global_log->info()
			<< "Forcing direct-pp communication scheme, as _forceDirectPP is set (probably by a child class)."
			<< std::endl;
		neighbourCommunicationScheme = "direct-pp";
	}

	std::string zonalMethod = "fs";
	std::string traversal = "c08"; // currently useless, as traversal is set elsewhere

	xmlconfig.changecurrentnode("../datastructure");
	xmlconfig.getNodeValue("traversalSelector", traversal);
	transform(traversal.begin(),
			  traversal.end(),
			  traversal.begin(),
			  ::tolower);
	// currently only checks, if traversal is valid - should check, if zonal method/traversal is valid
	if(traversal.find("hs") != std::string::npos || traversal.find("mp") != std::string::npos  || traversal.find("nt") != std::string::npos ) {
		zonalMethod = traversal;
	} else if(traversal.find("es") != std::string::npos){
		zonalMethod = "es";
	}
	else{
		Log::global_log->info() << "Defaulting to fs/c08" << std::endl;

		zonalMethod = "fs";
		traversal = "c08";
	}


	Log::global_log->info() << "variable zonalMethod is: " << zonalMethod << std::endl;
	// reset path
	xmlconfig.changecurrentnode(oldPath);

	// Specifies if the sequential fallback shall be used.
	bool useSequentialFallback = true;
	xmlconfig.getNodeValue("useSequentialFallback", useSequentialFallback);
	if (zonalMethod == "nt") {
		Log::global_log->info() << "Forcefully disabling sequential fallback, because Neutral Territory is used!" << std::endl;
		useSequentialFallback = false;
		Log::global_log->info() << "Enforcing direct-pp neighborcommunicationscheme, because NT is used!" << std::endl;
		neighbourCommunicationScheme = "direct-pp";
	}
	setCommunicationScheme(neighbourCommunicationScheme, zonalMethod);
	_neighbourCommunicationScheme->setSequentialFallback(useSequentialFallback);

	bool overlappingCollectives = false;
	xmlconfig.getNodeValue("overlappingCollectives", overlappingCollectives);
	if(overlappingCollectives) {
		Log::global_log->info() << "DomainDecompMPIBase: Using Overlapping Collectives" << std::endl;
#if MPI_VERSION >= 3
		_collCommunication = std::unique_ptr<CollectiveCommunicationInterface>(new CollectiveCommunicationNonBlocking());
#else
		Log::global_log->warning() << "DomainDecompMPIBase: Can not use overlapping collectives, as the MPI version is less than MPI 3." << std::endl;
#endif
		xmlconfig.getNodeValue("overlappingStartAtStep", _overlappingStartAtStep);
		Log::global_log->info() << "DomainDecompMPIBase: Overlapping Collectives start at step " << _overlappingStartAtStep
						   << std::endl;
	} else {
		Log::global_log->info() << "DomainDecompMPIBase: NOT Using Overlapping Collectives" << std::endl;
	}
}

int DomainDecompMPIBase::getNonBlockingStageCount() {
	return _neighbourCommunicationScheme->getCommDims();
}

void DomainDecompMPIBase::setCommunicationScheme(const std::string& scheme, const std::string& zonalMethod) {

	ZonalMethod* zonalMethodP = nullptr;

	// CommunicationScheme will delete the pointer
	if(zonalMethod=="fs") {
		zonalMethodP = new FullShell();
	} else if(zonalMethod=="es") {
		zonalMethodP = new EighthShell();
	} else if(zonalMethod=="hs") {
		zonalMethodP = new HalfShell();
	} else if(zonalMethod=="mp") {
		zonalMethodP = new Midpoint();
	} else if(zonalMethod=="nt") {
		zonalMethodP = new NeutralTerritory();
	} else {
		Log::global_log->error() << "DomainDecompMPIBase: invalid zonal method specified. Valid values are 'fs', 'es', 'hs', 'mp' and 'nt'"
				<< std::endl;
		Simulation::exit(1);
	}
	Log::global_log->info() << "Using zonal method: " << zonalMethod << std::endl;

	if (scheme=="direct") {
		Log::global_log->info() << "DomainDecompMPIBase: Using DirectCommunicationScheme without push-pull neighbors" << std::endl;
		_neighbourCommunicationScheme = std::make_unique<DirectNeighbourCommunicationScheme>(zonalMethodP, false);
	} else if(scheme=="direct-pp") {
		Log::global_log->info() << "DomainDecompMPIBase: Using DirectCommunicationScheme with push-pull neighbors" << std::endl;
		_neighbourCommunicationScheme = std::make_unique<DirectNeighbourCommunicationScheme>(zonalMethodP, true);
	} else if(scheme=="indirect") {
		Log::global_log->info() << "DomainDecompMPIBase: Using IndirectCommunicationScheme" << std::endl;
		_neighbourCommunicationScheme = std::make_unique<IndirectNeighbourCommunicationScheme>(zonalMethodP);
	} else {
		Log::global_log->error() << "DomainDecompMPIBase: invalid NeighbourCommunicationScheme specified. Valid values are 'direct' and 'indirect'"
				<< std::endl;
		Simulation::exit(1);
	}
}

unsigned DomainDecompMPIBase::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	std::vector<unsigned> moldistribution(_numProcs);
	MPI_CHECK(MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution.data(), 1, MPI_UNSIGNED, _comm));
	unsigned globalN = 0;
	for (int r = 0; r < _rank; r++)
		globalN += moldistribution[r];
	unsigned localNbottom = globalN;
	globalN += moldistribution[_rank];
	unsigned localNtop = globalN;
	for (int r = _rank + 1; r < _numProcs; r++)
		globalN += moldistribution[r];
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
				Log::global_log->error() << "IX is " << IX << " for rank 0, but " << recv << " for rank " << i << ".\n";
				MPI_Abort(_comm, 911);
			}
		}
		Log::global_log->info() << "IX = " << recv << " for all " << _numProcs << " ranks.\n";
	}
}

void DomainDecompMPIBase::assertDisjunctivity(ParticleContainer* moleculeContainer) const {

	// Put local IDs in vector
	unsigned long num_molecules = moleculeContainer->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
	std::vector<unsigned long> localParticleIds;
	localParticleIds.reserve(num_molecules);

	for (auto m = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); m.isValid(); ++m) {
		localParticleIds.push_back(m->getID());
	}

	// Check local duplicates
	std::unordered_set<long long> localDuplicates;
	std::unordered_set<long long> seenIds;

	for (const auto& pid : localParticleIds) {
		if (seenIds.find(pid) != seenIds.end()) {
			localDuplicates.insert(pid);
		} else {
			seenIds.insert(pid);
		}
	}

	// Distribute particle IDs to corresponding ranks
	std::vector<std::vector<long long>> sendBuffers(_numProcs);
	for (const auto& id : localParticleIds) {
		
		// Function to hash a particle ID to an MPI process rank
		int targetRank = id % _numProcs;

		sendBuffers[targetRank].push_back(id);
	}

	// Determine send counts and displacements
	std::vector<int> sendCounts(_numProcs, 0);
	std::vector<int> sendDispls(_numProcs, 0);
	for (int i = 0; i < _numProcs; ++i) {
		sendCounts[i] = sendBuffers[i].size();
	}

	// All-to-all communication: share the counts first
	std::vector<int> recvCounts(_numProcs);
	MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

	// Compute displacements for send and receive buffers
	int totalSend = 0, totalRecv = 0;
	std::vector<int> recvDispls(_numProcs, 0);
	for (int i = 0; i < _numProcs; ++i) {
		sendDispls[i] = totalSend;
		recvDispls[i] = totalRecv;
		totalSend += sendCounts[i];
		totalRecv += recvCounts[i];
	}

	// Prepare send and receive buffers
	std::vector<long long> sendBuffer(totalSend);
	std::vector<long long> recvBuffer(totalRecv);
	for (int i = 0; i < _numProcs; ++i) {
		std::copy(sendBuffers[i].begin(), sendBuffers[i].end(), sendBuffer.begin() + sendDispls[i]);
	}

	// All-to-all communication: exchange the particle IDs
	MPI_Alltoallv(sendBuffer.data(), sendCounts.data(), sendDispls.data(), MPI_LONG_LONG,
				  recvBuffer.data(), recvCounts.data(), recvDispls.data(), MPI_LONG_LONG, MPI_COMM_WORLD);

	// Check for duplicates in the received IDs
	std::unordered_set<long long> globalSeenIds;
	std::unordered_set<long long> globalDuplicates = localDuplicates; // Start with local duplicates

	for (const auto& pid : recvBuffer) {
		if (globalSeenIds.find(pid) != globalSeenIds.end()) {
			globalDuplicates.insert(pid);
		} else {
			globalSeenIds.insert(pid);
		}
	}

	// Collect and print global duplicates
	int globalDuplicateCount = globalDuplicates.size();
	std::vector<long long> globalDuplicateList(globalDuplicateCount);
	std::copy(globalDuplicates.begin(), globalDuplicates.end(), globalDuplicateList.begin());

	// Gather sizes of duplicate lists at root process
	std::vector<int> allDuplicateCounts(_numProcs);
	MPI_Gather(&globalDuplicateCount, 1, MPI_INT, allDuplicateCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Compute displacements for gathering duplicates at root
	std::vector<int> duplicateDispls(_numProcs, 0);
	int totalDuplicates = allDuplicateCounts[0];
	if (_rank == 0) {
		for (int i = 1; i < _numProcs; ++i) {
			duplicateDispls[i] = duplicateDispls[i - 1] + allDuplicateCounts[i - 1];
			totalDuplicates += allDuplicateCounts[i];
		}
	}

	// Gather all duplicates at root process
	std::vector<long long> allDuplicates(totalDuplicates);
	MPI_Gatherv(globalDuplicateList.data(), globalDuplicateCount, MPI_LONG_LONG,
				allDuplicates.data(), allDuplicateCounts.data(), duplicateDispls.data(),
				MPI_LONG_LONG, 0, MPI_COMM_WORLD);

	// Print global duplicates at root process
	if (_rank == 0) {
		if (totalDuplicates == 0) {
			Log::global_log->info() << "Data consistency checked: No duplicate IDs detected." << std::endl;
		} else {
			std::stringstream ss;
			for (const auto& pid : allDuplicates) {
				ss << pid << " ";
			}
			Log::global_log->error() << "Duplicate particle IDs found: " << ss.str() << std::endl;
			Log::global_log->error() << "Aborting because of duplicated particles." << std::endl;
			MPI_Abort(_comm, 1);
		}
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
		MessageType msgType, bool doHaloPositionCheck, bool removeRecvDuplicates) {

	Log::global_log->set_mpi_output_all();

	_neighbourCommunicationScheme->exchangeMoleculesMPI(moleculeContainer, domain, msgType, removeRecvDuplicates,
														this, doHaloPositionCheck);

	Log::global_log->set_mpi_output_root(0);
}



void DomainDecompMPIBase::exchangeForces(ParticleContainer* moleculeContainer, Domain* domain) {
	Log::global_log->set_mpi_output_all();

	// Using molecule exchange method with the force message type
	_neighbourCommunicationScheme->exchangeMoleculesMPI(moleculeContainer, domain, FORCES, false, this);

	Log::global_log->set_mpi_output_root(0);
}

size_t DomainDecompMPIBase::getTotalSize() {
	return DomainDecompBase::getTotalSize() + _neighbourCommunicationScheme->getDynamicSize()
			+ _collCommunication->getTotalSize();
}

void DomainDecompMPIBase::printDecomp(const std::string &filename, Domain *domain, ParticleContainer *particleContainer) {
	if (_rank == 0) {
		std::ofstream povcfgstrm(filename);
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " "
				   << domain->getGlobalLength(2) << std::endl;
		povcfgstrm << "rank boxMin_x boxMin_y boxMin_z boxMax_x boxMax_y boxMax_z Configuration" << std::endl;
		povcfgstrm.close();
	}

	std::stringstream localCellInfo;
	localCellInfo << _rank << " " << getBoundingBoxMin(0, domain) << " " << getBoundingBoxMin(1, domain) << " "
			<< getBoundingBoxMin(2, domain) << " " << getBoundingBoxMax(0, domain) << " "
			<< getBoundingBoxMax(1, domain) << " " << getBoundingBoxMax(2, domain) << " "
			<< particleContainer->getConfigurationAsString() << "\n";
	std::string localCellInfoStr = localCellInfo.str();

#ifdef ENABLE_MPI
	MPI_File fh;
	MPI_File_open(_comm, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_APPEND | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	uint64_t write_size = localCellInfoStr.size();
	uint64_t offset = 0;
	if (_rank == 0) {
		MPI_Offset file_end_pos;
		MPI_File_seek(fh, 0, MPI_SEEK_END);
		MPI_File_get_position(fh, &file_end_pos);
		write_size += file_end_pos;
		MPI_Exscan(&write_size, &offset, 1, MPI_UINT64_T, MPI_SUM, _comm);
		offset += file_end_pos;
	} else {
		MPI_Exscan(&write_size, &offset, 1, MPI_UINT64_T, MPI_SUM, _comm);
	}
	MPI_File_write_at(fh, static_cast<MPI_Offset>(offset), localCellInfoStr.c_str(), static_cast<int>(localCellInfoStr.size()), MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
#else
	std::ofstream povcfgstrm(filename.c_str(), std::ios::app);
	povcfgstrm << localCellInfoStr;
	povcfgstrm.close();
#endif
}

void DomainDecompMPIBase::printSubInfo(int offset) {
	std::stringstream offsetstream;
	for (int i = 0; i < offset; i++) {
		offsetstream << "\t";
	}
	Log::global_log->info() << offsetstream.str() << "own datastructures:\t" << sizeof(DomainDecompMPIBase) / 1.e6 << " MB" << std::endl;
	Log::global_log->info() << offsetstream.str() << "neighbourCommunicationScheme:\t\t" << _neighbourCommunicationScheme->getDynamicSize() / 1.e6 << " MB" << std::endl;
	Log::global_log->info() << offsetstream.str() << "collective Communication:\t\t" << _collCommunication->getTotalSize() / 1.e6 << " MB" << std::endl;

}

void DomainDecompMPIBase::printCommunicationPartners(std::string filename) const{
	_neighbourCommunicationScheme->printCommunicationPartners(filename);
}

void DomainDecompMPIBase::collCommAllreduceSumAllowPrevious() {
	if (global_simulation->getSimulationStep() >= _overlappingStartAtStep) {
		_collCommunication->allreduceSumAllowPrevious();
	} else {
		_collCommunication->allreduceSum();
	}
}
