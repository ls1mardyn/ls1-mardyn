/*
 * RedundancyResilience.h
 *
 *  Created on: 11 Jul 2018
 *      Author: tchipevn / Oliver Fernandes
 */

#ifdef ENABLE_MPI

#include "RedundancyResilience.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
using Log::global_log;

//temporary backup shit
#include <sstream>
using Log::global_log;

//pretty much default
ResilienceComm::ResilienceComm(int numProcs, int rank)
		: _numProcs(numProcs)
		, _rank(rank) {
}

//copy constructor
ResilienceComm::ResilienceComm(ResilienceComm const& rc)
		: _numProcs(rc._numProcs)
		, _rank(rc._rank) {
}

ResilienceComm::~ResilienceComm() {
	/*do nothing*/
}

int ResilienceComm::scatterBackupInfo(std::vector<int>& backupInfo, 
						  int const numberOfBackups,
						  size_t const sizePerRank,
	                      std::vector<int>& backing, 
	                      std::vector<int>& backedBy, 
	                      std::vector<int>& backingTags, 
	                      std::vector<int>& backedByTags) {
	size_t totalBytesRecv = sizePerRank*numberOfBackups*sizeof(int);
	std::vector<char> recvArray(totalBytesRecv);
	if (_rank == 0) {
		std::stringstream bkinf;
		// bkinf << "    RR: backupInfo: " << totalBytesRecv << "\n";
		// for (auto const& rnk : backupInfo) {
		// 	bkinf << rnk << ", ";
		// }
		// global_log->info() << bkinf.str() << std::endl;
		mardyn_assert(static_cast<uint>(sizePerRank*numberOfBackups*_numProcs) == backupInfo.size());
	}
	else {
		mardyn_assert(backupInfo.empty());
	}
	constexpr int const scatteringRank = 0;
	int mpi_error =	MPI_Scatter(reinterpret_cast<char*>(backupInfo.data()),
	                            totalBytesRecv,
	 			                MPI_CHAR,
				                recvArray.data(),
				                totalBytesRecv,
				                MPI_CHAR,
				                scatteringRank,
				                MPI_COMM_WORLD);
	mardyn_assert(mpi_error == MPI_SUCCESS);
	backing.resize(numberOfBackups);
	backedBy.resize(numberOfBackups);
	backingTags.resize(numberOfBackups);
	backedByTags.resize(numberOfBackups);
	auto backingAsChar = reinterpret_cast<char*>(backing.data());
	auto backedByAsChar = reinterpret_cast<char*>(backedBy.data());
	auto backingTagsAsChar = reinterpret_cast<char*>(backingTags.data());
	auto backedByTagsAsChar = reinterpret_cast<char*>(backedByTags.data());
	std::copy(recvArray.begin()+0*totalBytesRecv/sizePerRank, recvArray.begin()+1*totalBytesRecv/sizePerRank, backingAsChar);
	std::copy(recvArray.begin()+1*totalBytesRecv/sizePerRank, recvArray.begin()+2*totalBytesRecv/sizePerRank, backedByAsChar);
	std::copy(recvArray.begin()+2*totalBytesRecv/sizePerRank, recvArray.begin()+3*totalBytesRecv/sizePerRank, backingTagsAsChar);
	std::copy(recvArray.begin()+3*totalBytesRecv/sizePerRank, recvArray.begin()+4*totalBytesRecv/sizePerRank, backedByTagsAsChar);
	// global_log->info() << "    RR: Dumping scattered backup info: " << std::endl;
	// global_log->set_mpi_output_all();
	// std::stringstream bckd, bckBy, bckdTags, bckByTags;
	// for (int i=0; i<numberOfBackups; ++i) {
	// 	mardyn_assert(backing[i]<_numProcs);
	// 	mardyn_assert(backedBy[i]<_numProcs);
	// 	bckd << backing[i] << ", ";
	// 	bckBy << backedBy[i] << ", ";
	// 	bckdTags << backingTags[i] << ", ";
	// 	bckByTags << backedByTags[i] << ", ";
	// }
	// global_log->info() << "        Backed: " << bckd.str() << " Backed by: " << bckBy.str() << std::endl;
	// global_log->info() << "        Backed tags: " << bckdTags.str() << " Backed by tags: " << bckByTags.str() << std::endl;
	return 0;
}

int ResilienceComm::exchangeSnapshotSizes(
		std::vector<int>& backing,
		std::vector<int>& backedBy,
		std::vector<int>& backingTags,
		std::vector<int>& backedByTags,
		size_t const snapshotSize,
		std::vector<int>& backupDataSizes) {
	// send the size of this snapshot to all ranks backing it
	int src = -1;
	int dest = -1;
	int tag = -1;
	int status = MPI_ERR_UNKNOWN;
	for (size_t ib=0; ib<backedBy.size(); ++ib) {
		MPI_Request request=0;
		dest = backedBy[ib];
		tag = backedByTags[ib];
		// status = MPI_Isend(&snapshotSize, sizeof(snapshotSize), MPI_CHAR, dest, tag, MPI_COMM_WORLD, &request);
		status = MPI_Bsend(&snapshotSize, sizeof(snapshotSize), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
		mardyn_assert(status == MPI_SUCCESS);
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// setup the receiving buffers too for all ranks the current one is backing
	for (size_t ib=0; ib<backing.size(); ++ib) {
		MPI_Status recvStatus; 
		src = backing[ib];
		tag = backingTags[ib];
		void* target = &(backupDataSizes.data()[ib]);
		// status = MPI_Irecv(target, sizeof(snapshotSize), MPI_CHAR, src, tag, MPI_COMM_WORLD, &request);
		status = MPI_Recv(&(backupDataSizes.data()[ib]), sizeof(snapshotSize), MPI_CHAR, src, tag, MPI_COMM_WORLD, &recvStatus);
		mardyn_assert(status == MPI_SUCCESS);
	}
	mardyn_assert(status == MPI_SUCCESS);
	return 0;
}

int ResilienceComm::exchangeSnapshots(
		std::vector<int>& backing,
		std::vector<int>& backedBy,
		std::vector<int>& backingTags,
		std::vector<int>& backedByTags,
		std::vector<int>& backupDataSizes,
		std::vector<char>& sendData,
		std::vector<char>& recvData) {
	// send the snapshot to all ranks backing it
	int src = -1;
	int dest = -1;
	int tag = -1;
	int status = MPI_ERR_UNKNOWN;
	//prepare memory, create prefix sum
	std::vector<int> recvIndices(backupDataSizes.size());
	recvIndices[0] = 0; //first index always 0
	auto dstIt = recvIndices.begin()+1;
	auto srcIt = backupDataSizes.begin();
	while (dstIt != recvIndices.end()) {
		*dstIt = *(dstIt-1)+*srcIt;
		++dstIt; ++srcIt;
	}
	size_t const totalRecvSize = recvIndices.back()+*srcIt;
	recvData.resize(totalRecvSize);

	for (size_t ib=0; ib<backedBy.size(); ++ib) {
		dest = backedBy[ib];
		tag = backedByTags[ib];
		// global_log->info() << "    RR: Sending " << sendData.size()
		// 		<< " bytes to: " << dest << " using tag: " << tag << std::endl;
		status = MPI_Bsend(sendData.data(), sendData.size(), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
		mardyn_assert(status == MPI_SUCCESS);
	}
	// setup the receiving buffers too for all ranks the current one is backing
	for (size_t ib=0; ib<backing.size(); ++ib) {
		MPI_Status recvStatus;
		src = backing[ib];
		tag = backingTags[ib];
		size_t const recvIndex = recvIndices[ib];
		// global_log->info() << "    RR: Receiving " 
		// 		<< backupDataSizes[ib] << " bytes from " 
		// 		<< src << " at " 
		// 		<< recvIndices[ib] << " using tag: " << tag << std::endl;
		status = MPI_Recv(&recvData.data()[recvIndex], backupDataSizes[ib], MPI_CHAR, src, tag, MPI_COMM_WORLD, &recvStatus);
		mardyn_assert(status == MPI_SUCCESS);
		//verify a bunch of stuff
		int count;
		MPI_Get_count(&recvStatus, MPI_CHAR, &count);
		mardyn_assert(count == backupDataSizes[ib]);
		mardyn_assert(recvStatus.MPI_SOURCE == src);
		mardyn_assert(recvStatus.MPI_TAG == tag);
		mardyn_assert(status == MPI_SUCCESS);
	}
	return 0;
}

void RedundancyResilience::init(ParticleContainer* particleContainer,
								DomainDecompBase* domainDecomp,
								Domain* domain) {
	// instantiate the MPI communication object, and init it
	_comm = std::unique_ptr<ResilienceComm>(
			new ResilienceComm(domainDecomp->getNumProcs(),
	                domainDecomp->getRank()));
	// size per backup entry
	_sizePerRank = 4;
	// create a vector containing backup assignments (on root)
	std::vector< int > allBackupIds = _determineBackups(domainDecomp);
	// give each rank set of ids to backup (_backing) and tell each one which rank is backing it up (_backedBy)
	_comm->scatterBackupInfo(allBackupIds, 
			_numberOfBackups,
			_sizePerRank,
			_backing,
			_backedBy,
			_backingTags,
			_backedByTags);
}

void RedundancyResilience::readXML(XMLfileUnits& xmlconfig) {
	_backupInterval = 5;
	xmlconfig.getNodeValue("backupInterval", _backupInterval);
	global_log->info() << "    RR: Backup interval: " << _backupInterval << std::endl;

	_numberOfBackups = 1;
	xmlconfig.getNodeValue("numberOfBackups", _numberOfBackups);
	global_log->info() << "    RR: Number of ranks to back up: " << _numberOfBackups << std::endl;
}

///reset the particle container with a saved snapshot
void RedundancyResilience::beforeEventNewTimestep(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep) {
	// if (simstep != _restartAtIteration) {
	if (true) {
		return;
	}
	global_log->info() << "    RR: Resetting to " << _snapshot.getCurrentTime() << std::endl;
	Domain * domain = global_simulation->getDomain();

	// erase all current molecules
	particleContainer->clear();

	// fill new molecules
	particleContainer->addParticles(const_cast<std::vector<Molecule>& >(_snapshot.getMolecules()));

	// there should be no need to compute the forces again!
	// They are saved in the Molecule objects this time, due to the copy constructor.

	// Note that the forces, rotational moments are usually not saved in checkpoints and have to be recomputed in prepare_start()
	// so eventually, the following calls may be necessary:
	// * Simulation::updateParticleContainerAndDecomposition,
	// * ParticleContainer::traverseCells(cellProcessor)
	// * Simulation::updateForces()

	// set globals
	global_simulation->setSimulationTime(_snapshot.getCurrentTime());
	domain->setGlobalTemperature(_snapshot.getTemperature());

	mardyn_assert(_snapshot.getGlobalNumberOfMolecules() == domain->getglobalNumMolecules());
	mardyn_assert(domainDecomp->getRank() == _snapshot.getRank());
}

void RedundancyResilience::endStep(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	if (simstep % _backupInterval != 0) {
		return;
	}
	global_log->info() << "    RR: Creating in-memory backups... " << std::endl;
	_saveLocalSnapshot(particleContainer, domainDecomp, domain, simstep);
	std::vector<char> snapshotDataAsBytes = _serializeSnapshot();
	std::vector<int> backupDataSizes(_numberOfBackups);
	std::vector<char> backupData;
	_comm->exchangeSnapshotSizes(
			_backing,
			_backedBy,
			_backingTags,
			_backedByTags,
			snapshotDataAsBytes.size(),
			backupDataSizes);
	_comm->exchangeSnapshots(
			_backing,
			_backedBy,
			_backingTags,
			_backedByTags,
			backupDataSizes,
			snapshotDataAsBytes,
			backupData);
	global_log->info() << "    RR: Snapshots exchanged for timestep " << simstep << std::endl;
	_storeSnapshots(backupDataSizes, backupData);
}

void RedundancyResilience::_saveLocalSnapshot(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	// put the molecules in the buffer
	_snapshot.clearMolecules();
	for (auto m = particleContainer->iterator(); m.isValid(); ++m) {
		_snapshot.addMolecule(*m);
	}
	//set time, global number of molecules and temperature
	_snapshot.setCurrentTime(global_simulation->getSimulationTime());
	_snapshot.setGlobalNumberOfMolecules(domain->getglobalNumMolecules());
	_snapshot.setTemperature(domain->getGlobalCurrentTemperature());
	_snapshot.setRank(domainDecomp->getRank());
}

std::vector<char> RedundancyResilience::_serializeSnapshot(void) const {
	auto const rank = _snapshot.getRank();
	auto const currentTime = _snapshot.getCurrentTime();

	std::vector<char> byteData;
	byteData.insert(byteData.end(), 
			reinterpret_cast<char const*>(&rank), 
			reinterpret_cast<char const*>(&rank)+sizeof(rank));
	mardyn_assert(byteData.size() == sizeof(rank));
	byteData.insert(byteData.end(), 
	        reinterpret_cast<char const*>(&currentTime), 
			reinterpret_cast<char const*>(&currentTime)+sizeof(currentTime));
	mardyn_assert(byteData.size() == sizeof(rank)+sizeof(currentTime));
	//append 0,...,rank as char to generate different sized data for each rank while encoding some info
#warning serializing and deserializing is generating fake data
	for (char fakeData = 0; fakeData<rank+1; ++fakeData) {
		byteData.push_back(fakeData);
	}
	mardyn_assert(byteData.size() == sizeof(rank)+sizeof(currentTime)+static_cast<size_t>(rank)+1);
	return byteData;
}

void RedundancyResilience::_storeSnapshots(std::vector<int>& backupDataSizes, std::vector<char>& backupData) {
	auto snapshotStartIt = backupData.begin();
	auto snapshotSizeIt = backupDataSizes.begin();
	_backupSnapshots.clear();
	while (snapshotStartIt != backupData.end()) {
		Snapshot newSnapshot;
		assert(snapshotSizeIt != backupDataSizes.end());
		auto snapshotEndIt = snapshotStartIt+*snapshotSizeIt;
		snapshotStartIt = _deserializeSnapshot(snapshotStartIt, snapshotEndIt, newSnapshot);
		++snapshotSizeIt;
		_backupSnapshots.push_back(newSnapshot);
	}
	mardyn_assert(_backupSnapshots.size() == static_cast<size_t>(_numberOfBackups));
}

std::vector<char>::iterator RedundancyResilience::_deserializeSnapshot(std::vector<char>::iterator const snapshotStart,
        std::vector<char>::iterator const snapshotEnd, Snapshot& newSnapshot) {
	int rank;
	double currentTime;
	auto valueStart = snapshotStart;
	auto valueEnd = snapshotStart+sizeof(rank);
	std::copy(valueStart, valueEnd, reinterpret_cast<char*>(&rank));
	newSnapshot.setRank(rank);
	valueStart = valueEnd;
	valueEnd = valueStart+sizeof(currentTime);
	std::copy(valueStart, valueEnd, reinterpret_cast<char*>(&currentTime));
	valueStart = valueEnd;
	valueEnd = snapshotEnd;
	// deserialize the fake data, used for debug purposes
	std::vector<char> fakeData(valueEnd - valueStart);
	std::copy(valueStart, valueEnd, fakeData.begin());
	_validateFakeData(rank, fakeData);
	mardyn_assert(valueEnd == snapshotEnd);
	return snapshotEnd;
}

bool RedundancyResilience::_validateFakeData(int const rank, std::vector<char>& fakeData) {
	mardyn_assert(fakeData.size() == static_cast<size_t>(rank+1));
	mardyn_assert(*fakeData.begin() == 0);
	mardyn_assert(*(fakeData.end()-1) == rank);
	int i=0;
	for (auto fDit=fakeData.begin(); fDit!=fakeData.end(); ++fDit, ++i) {
		mardyn_assert(*fDit == i);
	}
	return true;
}

std::vector<int> RedundancyResilience::_determineBackups(DomainDecompBase const* domainDecomp) {
	std::vector<int> backupInfo;
	backupInfo.clear();
	if (domainDecomp->getRank() == 0) {
		global_log->info() << "    RR: Determining new backup rank distribution. " << std::endl;
		int const numRanks = domainDecomp->getNumProcs();
		backupInfo.resize(_sizePerRank*numRanks*_numberOfBackups);
		// generate the pairs rank->backed up rank. Also fill out the reverse.
		// and enumerate the association the tag for both
		int tag = 0;
		for (int ib=0; ib<_numberOfBackups; ++ib) {
			for (int rank=0; rank<numRanks; ++rank) {
				// determine which node backs up what by some clever scheme 
				// (here: stupid scheme which simply selects ranks {n+i | i=1,...,numberOfBackups} for backing rank n)
				int newBackup = rank+ib+1;
				if (newBackup >= numRanks) newBackup=newBackup%numRanks;
				mardyn_assert(newBackup < numRanks);
				// => newBackup will be backed up by rank
				int const backupIdx = rank*_numberOfBackups*_sizePerRank + 0*_numberOfBackups + ib;
				int const backedByIdx = newBackup*_numberOfBackups*_sizePerRank + 1*_numberOfBackups + ib;
				int const tagBackupIdx = backupIdx + 2*_numberOfBackups;
				int const tagBackedByIdx = backedByIdx + 2*_numberOfBackups;
				backupInfo[backupIdx] = newBackup;
				backupInfo[backedByIdx] = rank;
				backupInfo[tagBackupIdx] = tag;
				backupInfo[tagBackedByIdx] = tag;
				++tag;
			}
		}
		mardyn_assert(tag == numRanks*_numberOfBackups);
		mardyn_assert(backupInfo.size() == static_cast<unsigned int>(_sizePerRank*numRanks*_numberOfBackups));
	}

	// std::stringstream bkinf;
	// bkinf << "    RR: backupInfo: \n";
	// for (auto const& rnk : backupInfo) {
	// 	bkinf << rnk << ", ";
	// }
	// global_log->info() << bkinf.str() << std::endl;

	return backupInfo;
}
#endif /*ENABLE_MPI */