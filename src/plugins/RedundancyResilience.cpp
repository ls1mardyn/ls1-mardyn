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
	_comm = std::unique_ptr<ResilienceComm>(new ResilienceComm(domainDecomp->getNumProcs(),
	                                                           domainDecomp->getRank()));
	// create a vector containing backup assignments (on root)
	std::vector< int > allBackupIds = _determineBackups(domainDecomp);
	// give each rank set of ids to backup (_backing) and tell each one which rank is backing it up (_backedBy)
	_comm->scatterBackupInfo(allBackupIds, _backing, _backedBy, _numberOfBackups);
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
	global_log->info() << "    RR: Creating in-memory backups " << simstep << std::endl;
	_saveLocalSnapshot(particleContainer, domainDecomp, domain, simstep);
	std::vector<char> snapshotDataAsBytes = _serializeSnapshot();
	std::vector<int> backupDataSizes(_numberOfBackups);
	std::vector<char> backupData;
	_comm->exchangeSnapshotSizes(_backing,
					    	     _backedBy,
								 snapshotDataAsBytes.size(),
						         backupDataSizes);
	_comm->exchangeSnapshots(_backing,
					    	 _backedBy,
						     backupDataSizes,
						     snapshotDataAsBytes,
						     backupData);
	_storeSnapshots(backupDataSizes, backupData);	
}

void RedundancyResilience::_saveLocalSnapshot(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) {
	// put the molecules in the buffer
	_snapshot.clearMolecules();
	for (ParticleIterator m = particleContainer->iterator(); m.hasNext(); m.next()) {
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

	size_t totalBytes = sizeof(rank)+sizeof(currentTime);
	std::vector<char> byteData;
	byteData.insert(byteData.end(), 
	                reinterpret_cast<char const*>(&rank), 
					reinterpret_cast<char const*>(&rank)+sizeof(rank));
	byteData.insert(byteData.end(), 
	                reinterpret_cast<char const*>(&currentTime), 
					reinterpret_cast<char const*>(&currentTime)+sizeof(currentTime));
	mardyn_assert(byteData.size() == totalBytes);
	return byteData;
}

void RedundancyResilience::_storeSnapshots(std::vector<int>& backupDataSizes, std::vector<char>& backupData) {
	auto snapshotStart = backupData.begin();
	auto snapshotSize = backupDataSizes.begin();
	_backupSnapshots.clear();
	while (snapshotStart != backupData.end()) {
		Snapshot newSnapshot;
		assert(snapshotSize != backupDataSizes.end());
		auto snapshotEnd = snapshotStart+*snapshotSize;
		snapshotStart = _deserializeSnapshot(snapshotStart, snapshotEnd, newSnapshot);
		++snapshotSize;
		_backupSnapshots.push_back(newSnapshot);
	}
	global_log->info() << "_backupSnapshots.size(): " << _backupSnapshots.size() << " numberOfBackups: " <<  static_cast<size_t>(_numberOfBackups) << std::endl;
	mardyn_assert(_backupSnapshots.size() == static_cast<size_t>(_numberOfBackups));
}

std::vector<char>::iterator RedundancyResilience::_deserializeSnapshot(std::vector<char>::iterator snapshotStart,
        std::vector<char>::iterator snapshotEnd, Snapshot& snapshot) {
	int rank;
	double currentTime;
	auto valueStart = snapshotStart;
	auto valueEnd = snapshotStart+sizeof(rank);
	std::copy(valueStart, valueEnd, reinterpret_cast<char*>(&rank));
	valueStart = valueEnd;
	valueEnd = valueStart+sizeof(currentTime);
	std::copy(snapshotStart, snapshotStart+sizeof(rank), reinterpret_cast<char*>(&rank));
	mardyn_assert(valueEnd == snapshotEnd);
	return snapshotEnd;
}

std::vector<int> RedundancyResilience::_determineBackups(DomainDecompBase const* domainDecomp) {
	std::vector<int> backupInfo;
	backupInfo.clear();
	if (domainDecomp->getRank() == 0) {
		global_log->info() << "    RR: Determining new backup rank distribution. " << std::endl;
		int const numRanks = domainDecomp->getNumProcs();
		backupInfo.resize(2*numRanks*_numberOfBackups);
		// ib denotes the i-th backup
		for (int ib=0; ib<_numberOfBackups; ++ib) {
			for (int rank=0; rank<numRanks; ++rank) {
				// determine which node backs up what by some clever scheme 
				// (here: stupid scheme which simply selects ranks {n+i | i=1,...,numberOfBackups} for rank n)
				int newBackup = rank+ib+1;
				if (newBackup >= numRanks) newBackup=newBackup%numRanks;
				mardyn_assert(newBackup < numRanks);
				// => newBackup will be backed up by rank
				int const backupIdx = rank*_numberOfBackups*2+ib;
				int const backedByIdx = newBackup*_numberOfBackups*2+_numberOfBackups+ib;
				backupInfo[backupIdx] = newBackup;
				backupInfo[backedByIdx] = rank;
			}
		}
		mardyn_assert(backupInfo.size() == static_cast<unsigned int>(2*numRanks*_numberOfBackups));
	}
	return backupInfo;
}
#endif /*ENABLE_MPI */