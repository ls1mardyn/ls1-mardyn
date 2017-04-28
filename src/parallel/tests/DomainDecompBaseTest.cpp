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

void DomainDecompBaseTest::testExchangeMoleculesSimple() {
	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,1,1,1,0,false);
	components.push_back(dummyComponent);

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
	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,1,1,1,0,false);
	components.push_back(dummyComponent);

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
