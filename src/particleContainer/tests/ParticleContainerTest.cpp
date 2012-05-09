/*
 * ParticleContainerTest.cpp
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#include "ParticleContainerTest.h"
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
	Molecule dummyMolecule1(1, 0, 1.0,1.0,1.0,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule1);
	Molecule dummyMolecule2(2, 0, 2.0,2.0,2.0,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule2);
	Molecule dummyMolecule3(3, 0, 3.0,3.0,3.0,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule3);
	Molecule dummyMolecule4(4, 0, 5.1,5.1,5.1,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule4);
}

void ParticleContainerTest::testInsertion(ParticleContainer* container) {
	setupMolecules(container);
	ASSERT_EQUAL(4ul, container->getNumberOfParticles());
	Molecule dummyMolecule5(0, 0, 7.7,7.1,7.1,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule5);
	Molecule dummyMolecule6(0, 0, 11.1,11.1,11.1,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule6);
	ASSERT_EQUAL(6ul, container->getNumberOfParticles());
}


void ParticleContainerTest::testMoleculeIteration(ParticleContainer* container) {
	setupMolecules(container);
	Molecule* molecule = container->begin();
	int moleculeCount = 0;
	while (molecule != container->end()) {
		molecule = container->next();
		moleculeCount++;
	}
	ASSERT_EQUAL(4, moleculeCount);
	ASSERT_EQUAL(4ul, container->getNumberOfParticles());

	Molecule dummyMolecule5(5, 0, 7.7,7.1,7.1,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule5);
	Molecule dummyMolecule6(6, 0, 11.1,11.1,11.1,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule6);

	moleculeCount = 0;
	bool ids[] = {false, false, false, false, false, false};

	molecule = container->begin();
	while (molecule != container->end()) {
		ids[molecule->id() - 1] = true;
		test_log->debug() << "Visited Molecule with id " << molecule->id() << std::endl;
		molecule = container->next();
		moleculeCount++;
	}
	ASSERT_EQUAL(6, moleculeCount);
	ASSERT_EQUAL(6ul, container->getNumberOfParticles());
	for (int i = 0; i < 6; i++) {
		ASSERT_TRUE(ids[i]);
	}
}


void ParticleContainerTest::testUpdateAndDeleteOuterParticles(ParticleContainer* container) {
	setupMolecules(container);

	// iterate over molecules and move
	Molecule* molecule = container->begin();
	int moleculeCount = 0;
	while (molecule != container->end()) {
		moleculeCount++;

		if (molecule->id() == 1) {
			molecule->setr(0, -0.2);
		} else if (molecule->id() == 3) {
			molecule->setr(1, 9.0);
			molecule->setr(2, 9.4);
		}
		molecule = container->next();
	}

	ASSERT_EQUAL(4, moleculeCount);
	ASSERT_EQUAL(4ul, container->getNumberOfParticles());

	// add molecules and call update
	Molecule dummyMolecule5(5, 0, 7.7,7.1,7.1,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule5);
	Molecule dummyMolecule6(6, 0, 11.1,11.1,11.1,0,0,0, 0, 0, 0, 0, 0, 0, 0, &_components);
	container->addParticle(dummyMolecule6);
	container->update();

	moleculeCount = 0;
	bool ids[] = {false, false, false, false, false, false};

	molecule = container->begin();
	while (molecule != container->end()) {
		ids[molecule->id() - 1] = true;
		test_log->debug() << "Visited Molecule with id " << molecule->id() << std::endl;
		molecule = container->next();
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

	molecule = container->begin();
	moleculeCount = 0;
	while (molecule != container->end()) {
		ids[molecule->id() - 1] = true;
		test_log->debug() << "Visited Molecule with id " << molecule->id() << std::endl;
		molecule = container->next();
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
