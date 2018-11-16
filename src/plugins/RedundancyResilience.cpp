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
#pragma message "serializing and deserializing is generating fake data"
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