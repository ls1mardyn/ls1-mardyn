#include "ParticleContainerTest.h"

#include <set>

#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"

ParticleContainerTest::ParticleContainerTest() {
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,1,1,1,0,false);
	_components.push_back(dummyComponent);
}

ParticleContainerTest::~ParticleContainerTest() {
}

void ParticleContainerTest::setupMolecules(ParticleContainer* container) {
	//(id, cid, x, y, z, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz)
	Molecule dummyMolecule1(1, &_components[0], 1.0,1.0,1.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	container->addParticle(dummyMolecule1);
	Molecule dummyMolecule2(2, &_components[0], 2.0,2.0,2.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	container->addParticle(dummyMolecule2);
	Molecule dummyMolecule3(3, &_components[0], 3.0,3.0,3.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	container->addParticle(dummyMolecule3);
	Molecule dummyMolecule4(4, &_components[0], 5.1,5.1,5.1,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	container->addParticle(dummyMolecule4);
}

void ParticleContainerTest::testInsertion(ParticleContainer* container) {
	setupMolecules(container);
	ASSERT_EQUAL(4ul, container->getNumberOfParticles());
	Molecule dummyMolecule5(0, &_components[0], 7.7,7.1,7.1,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	ASSERT_EQUAL(true, container->addParticle(dummyMolecule5));
	Molecule dummyMolecule6(0, &_components[0], 11.1,11.1,11.1,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	ASSERT_EQUAL(true, container->addParticle(dummyMolecule6));
	ASSERT_EQUAL(6ul, container->getNumberOfParticles());
}


void ParticleContainerTest::testMoleculeIteration(ParticleContainer* container) {
	setupMolecules(container);
	unsigned long moleculeCount = 0;
	std::set<unsigned long> ids;
	for(auto moleculeIter = container->iteratorBegin(); moleculeIter != container->iteratorEnd(); ++moleculeIter) {
		test_log->debug() << "Visited Molecule with id " << moleculeIter->id() << std::endl;
		ids.insert(moleculeIter->id());
		moleculeCount++;
	}
	ASSERT_EQUAL(container->getNumberOfParticles(), moleculeCount); // check if iterator catches all molecules
	ASSERT_EQUAL(moleculeCount, ids.size()); // check for duplicate molecules
}


void ParticleContainerTest::testUpdateAndDeleteOuterParticles(ParticleContainer* container) {
	setupMolecules(container);

	// iterate over molecules and move
	ParticleIterator molecule = container->iteratorBegin();
	int moleculeCount = 0;
	while (molecule != container->iteratorEnd()) {
		moleculeCount++;

		if (molecule->id() == 1) {
			molecule->setr(0, -0.2);
		} else if (molecule->id() == 3) {
			molecule->setr(1, 1.0);
			molecule->setr(2, 2.4);
		}
		++molecule;
	}

	ASSERT_EQUAL(4, moleculeCount);
	ASSERT_EQUAL(4ul, container->getNumberOfParticles());

	// add molecules and call update
	Molecule dummyMolecule5(5, &_components[0], 7.7,7.1,7.1,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	container->addParticle(dummyMolecule5);
	Molecule dummyMolecule6(6, &_components[0], 11.1,11.1,11.1,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	container->addParticle(dummyMolecule6);
	container->update();

	moleculeCount = 0;
	bool ids[] = {false, false, false, false, false, false};

	molecule = container->iteratorBegin();
	while (molecule != container->iteratorEnd()) {
		ids[molecule->id() - 1] = true;
		test_log->debug() << "Visited Molecule with id " << molecule->id() << std::endl;
		++molecule;
		moleculeCount++;
	}
	ASSERT_EQUAL(6, moleculeCount);
	ASSERT_EQUAL(6ul, container->getNumberOfParticles());
	for (int i = 0; i < 6; i++) {
		ASSERT_TRUE(ids[i]);
	}

	// delete outer particles (i.e. particles 1 and 6) and check
	container->deleteOuterParticles();
	ASSERT_EQUAL(4ul, container->getNumberOfParticles());
	for (int i = 0; i < 6; i++) ids[i] = false;

	molecule = container->iteratorBegin();
	moleculeCount = 0;
	while (molecule != container->iteratorEnd()) {
		ids[molecule->id() - 1] = true;
		test_log->debug() << "Visited Molecule with id " << molecule->id() << std::endl;
		++molecule;
		moleculeCount++;
	}
	ASSERT_EQUAL(4, moleculeCount);

	for (int i = 0; i < 6; i++) {
		if (i == 0 || i == 5) {
			ASSERT_TRUE(!ids[i]);
		} else {
			ASSERT_TRUE(ids[i]);
		}
	}
}
