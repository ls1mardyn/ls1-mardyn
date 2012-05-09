/*
 * DomainDecompBaseTest.cpp
 *
 * @Date: 09.05.2012
 * @Author: eckhardw
 */

#include "parallel/tests/DomainDecompBaseTest.h"
#include "parallel/DomainDecompDummy.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Component.h"
#include "molecules/Molecule.h"

TEST_SUITE_REGISTRATION(DomainDecompBaseTest);

DomainDecompBaseTest::DomainDecompBaseTest() {
}

DomainDecompBaseTest::~DomainDecompBaseTest() {
}

void DomainDecompBaseTest::testExchangeMolecules() {
	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,1,1,1,0,false);
	components.push_back(dummyComponent);

	// make sure we have a DomainDecompDummy
	delete _domainDecomposition;
	_domainDecomposition = new DomainDecompDummy();
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "DomainDecompBase.inp", 17.0);

	unsigned int count[3] = {0};
	ASSERT_EQUAL(3ul, container->getNumberOfParticles());

	Molecule* m = container->begin();
	while ( m != container->end()) {
		count[m->id()]++;
		m = container->next();
	}

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(1u, count[i]);
		count[i] = 0;
	}

	// after the exchange, there have to be 21 copies in the halos, i.e. 24 molecules in total
	_domainDecomposition->exchangeMolecules(container, components, _domain);
	ASSERT_EQUAL(24ul, container->getNumberOfParticles());

	m = container->begin();
	while(m != container->end()) {
		count[m->id()]++;
		m = container->next();
	}

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(8u, count[i]);
	}
}
