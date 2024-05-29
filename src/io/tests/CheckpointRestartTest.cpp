/*
 * CheckpointRestartTest.cpp
 *
 * Check whether a checkpoint can be successfully read again.
 *
 *  Created on: 11.08.2016
 *      Author: seckler
 */

#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include <iostream>

#include "io/tests/CheckpointRestartTest.h"


#if !defined(ENABLE_REDUCED_MEMORY_MODE)
TEST_SUITE_REGISTRATION(CheckpointRestartTest);
#else
#pragma message "Compilation Info: CheckpointRestartTest disabled in reduced memory mode."
#endif

/*
 * This tests if a written checkpoint can successfully be read again using the ascii writer.
 */
void CheckpointRestartTest::testCheckpointRestartASCII() {
	testCheckpointRestart(false);
}

/*
 * This tests if a written checkpoint can successfully be read again using binary.
 */
void CheckpointRestartTest::testCheckpointRestartBinary() {
	testCheckpointRestart(true);
}

/*
 * Actual test if a written checkpoint can successfully be read again.
 */
void CheckpointRestartTest::testCheckpointRestart(bool binary) {
	constexpr double cutoff = 10.5;
	ParticleContainer* particleContainer
		= initializeFromFile(ParticleContainerFactory::LinkedCell, "VectorizationMultiComponentMultiPotentials_50_molecules.inp", cutoff);
	auto initialParticleCount = getGlobalParticleNumber(particleContainer);

	std::string filename = binary ? "restart.test" : "restart.test.dat";
	_domain->writeCheckpoint(getTestDataFilename(filename, false), particleContainer, _domainDecomposition,
	                         0., binary);

	delete particleContainer;

	ParticleContainer* particleContainer2
		= initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff, binary);

	auto restartedParticleCount = getGlobalParticleNumber(particleContainer2);
	ASSERT_EQUAL(initialParticleCount, restartedParticleCount);
	delete particleContainer2;
}

unsigned long CheckpointRestartTest::getGlobalParticleNumber(ParticleContainer* particleContainer){
	unsigned long localParticleCount = particleContainer->getNumberOfParticles();
#ifdef ENABLE_PERSISTENT
	auto collComm = make_CollCommObj_AllreduceAdd(_domainDecomposition->getCommunicator(), localParticleCount);
	collComm.persistent();
	unsigned long globalNumParticles;
	collComm.get(globalNumParticles);
#else
	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendUnsLong(localParticleCount);
	_domainDecomposition->collCommAllreduceSum();
	auto globalNumParticles = _domainDecomposition->collCommGetUnsLong();
	_domainDecomposition->collCommFinalize();
#endif
	return globalNumParticles;
}
