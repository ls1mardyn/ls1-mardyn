/*
 * DomainDecompositionTest.cpp
 *
 * @Date: 10.05.2012
 * @Author: eckhardw
 */

#include "DomainDecompositionTest.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecomposition.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Component.h"
#include "molecules/Molecule.h"
#include "Domain.h"

TEST_SUITE_REGISTRATION(DomainDecompositionTest);

DomainDecompositionTest::DomainDecompositionTest() {
}

DomainDecompositionTest::~DomainDecompositionTest() {
}

void DomainDecompositionTest::testNoDuplicatedParticlesFilename(const char * filename, double cutoff) {
	// original pointer will be deleted by tearDown()
	_domainDecomposition = new DomainDecomposition();

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);
	int numMols = container->getNumberOfParticles();

	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendInt(numMols);
	_domainDecomposition->collCommAllreduceSum();
	numMols = _domainDecomposition->collCommGetInt();
	_domainDecomposition->collCommFinalize();

	_domainDecomposition->balanceAndExchange(0.,true, container, _domain);
	container->deleteOuterParticles();

	int newNumMols = container->getNumberOfParticles();

	_domainDecomposition->collCommInit(1);
	_domainDecomposition->collCommAppendInt(newNumMols);
	_domainDecomposition->collCommAllreduceSum();
	newNumMols = _domainDecomposition->collCommGetInt();
	_domainDecomposition->collCommFinalize();
	std::cout << "old: " << numMols << " molecules; "<< " new: " << newNumMols << " molecules "<< std::endl;
	ASSERT_EQUAL(numMols, newNumMols);

	delete _domainDecomposition;
}

void DomainDecompositionTest::testNoDuplicatedParticles() {
	testNoDuplicatedParticlesFilename("H20_NaBr_0.01_T_293.15_DD.inp", 5.0);
}


void DomainDecompositionTest::testExchangeMolecules1Proc() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}

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
}
