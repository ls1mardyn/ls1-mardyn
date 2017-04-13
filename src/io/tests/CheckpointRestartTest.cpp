/*
 * inputFileTest.cpp
 *
 * There are two problems
 * - How to derive components.
 * - error with the use of molecule
 *
 *  Created on: 01.05.2012
 *      Author: yutaka
 */

#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include <iostream>
#include <sstream>

#include "io/tests/CheckpointRestartTest.h"

using namespace std;

class MDGenerator;

#ifndef MARDYN_WR
TEST_SUITE_REGISTRATION(CheckpointRestartTest);
#else
#warning "CheckpointRestartTest disabled in MARDYN_WR mode."
#endif


CheckpointRestartTest::CheckpointRestartTest() {
}

CheckpointRestartTest::~CheckpointRestartTest() {
}

/*
 * testRemoveMomentum tests if removeMomentum in MDGenerator works properly or not.
 */
void CheckpointRestartTest::testCheckpointRestart() {
	ParticleContainer* particleContainer
		= initializeFromFile(ParticleContainerFactory::LinkedCell, "VectorizationMultiComponentMultiPotentials_50_molecules.inp", 1.5);

	_domain->writeCheckpoint(getTestDataFilename("restart.test.xdr", false), particleContainer, _domainDecomposition, 0.);

	delete particleContainer;

	ParticleContainer* particleContainer2
			= initializeFromFile(ParticleContainerFactory::LinkedCell, "restart.test.xdr", 1.5);

	delete particleContainer2;
}
