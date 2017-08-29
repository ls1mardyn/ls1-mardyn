/*
 * DomainDecompBaseTest.cpp
 *
 * @Date: 09.05.2012
 * @Author: eckhardw
 */

#include "parallel/tests/DomainDecompBaseTest.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Component.h"
#include "molecules/Molecule.h"

TEST_SUITE_REGISTRATION(DomainDecompBaseTest);

DomainDecompBaseTest::DomainDecompBaseTest() {
}

DomainDecompBaseTest::~DomainDecompBaseTest() {
}

void DomainDecompBaseTest::testNoDuplicatedParticlesFilename(const char * filename, double cutoff) {
	// original pointer will be deleted by tearDown()
	_domainDecomposition = new DomainDecompBase();

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);
	int numMols = container->getNumberOfParticles();

	_domainDecomposition->exchangeMolecules(container, _domain);
	container->deleteOuterParticles();

	int newNumMols = container->getNumberOfParticles();
//	_domain->writeCheckpoint("dump.txt", container, _domainDecomposition);
	ASSERT_EQUAL(numMols, newNumMols);

	delete _domainDecomposition;
}

void DomainDecompBaseTest::testNoDuplicatedParticles() {
	testNoDuplicatedParticlesFilename("H20_NaBr_0.01_T_293.15.inp", 27.0);
}

void DomainDecompBaseTest::testExchangeMoleculesSimple() {

	// make sure we have a DomainDecompBase
	_domainDecomposition = new DomainDecompBase();
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "LinkedCells.inp", 1.0);

	unsigned int count[8] = {0};
	ASSERT_EQUAL(8ul, container->getNumberOfParticles());

	ParticleIterator m = container->iteratorBegin();
	while ( m != container->iteratorEnd()) {
		count[m->id()]++;
		++m;
	}

	for (int i = 0; i < 8; i++) {
		ASSERT_EQUAL(1u, count[i]);
		count[i] = 0;
	}

	// after the exchange, every molecule should be replicated 7 times, i.e. there have to be 64 molecules in total
	_domainDecomposition->exchangeMolecules(container, _domain);
	ASSERT_EQUAL(64ul, container->getNumberOfParticles());

	m = container->iteratorBegin();
	while(m != container->iteratorEnd()) {
		count[m->id()]++;
		++m;
	}

	for (int i = 0; i < 8; i++) {
		ASSERT_EQUAL(8u, count[i]);
	}

	delete container;
	delete _domainDecomposition;
}

void DomainDecompBaseTest::testExchangeMolecules() {

	// make sure we have a DomainDecompBase
	_domainDecomposition = new DomainDecompBase();
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "DomainDecompBase.inp", 17.0);

	unsigned int count[3] = {0};
	ASSERT_EQUAL(3ul, container->getNumberOfParticles());

	ParticleIterator m = container->iteratorBegin();
	while ( m != container->iteratorEnd()) {
		count[m->id()]++;
		++m;
	}

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(1u, count[i]);
		count[i] = 0;
	}

	// after the exchange, there have to be 21 copies in the halos, i.e. 24 molecules in total
	_domainDecomposition->exchangeMolecules(container, _domain);
	ASSERT_EQUAL(24ul, container->getNumberOfParticles());

	m = container->iteratorBegin();
	while(m != container->iteratorEnd()) {
		count[m->id()]++;
		++m;
	}

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(8u, count[i]);
	}

	delete container;
	delete _domainDecomposition;
}
