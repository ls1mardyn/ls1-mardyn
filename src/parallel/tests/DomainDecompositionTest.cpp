/*
 * DomainDecompositionTest.cpp
 *
 * @Date: 10.05.2012
 * @Author: eckhardw
 */

#include "DomainDecompositionTest.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Component.h"
#include "molecules/Molecule.h"

TEST_SUITE_REGISTRATION(DomainDecompositionTest);

DomainDecompositionTest::DomainDecompositionTest() {
}

DomainDecompositionTest::~DomainDecompositionTest() {
}


void DomainDecompositionTest::testExchangeMolecules1Proc() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}

	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,1,1,1,0,false);
	components.push_back(dummyComponent);

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
	_domainDecomposition->exchangeMolecules(container, _domain);
	ASSERT_EQUAL(24ul, container->getNumberOfParticles());

	m = container->begin();
	while(m != container->end()) {
		count[m->id()]++;
		m = container->next();
	}

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(8u, count[i]);
	}

	delete container;
}
