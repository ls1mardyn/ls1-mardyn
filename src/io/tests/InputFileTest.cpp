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
#include "molecules/Molecule.h"
#include "molecules/Component.h"
#include <iostream>
#include <sstream>

#include "io/tests/InputFileTest.h"
#include "../tools/gui/generators/MDGenerator.h"

using namespace std;
#ifdef SUPPORT_GENERATOR

class MDGenerator;

TEST_SUITE_REGISTRATION(InputFileTest);

InputFileTest::InputFileTest() {
}

InputFileTest::~InputFileTest() {
}

/*
 * testRemoveMomentum tests if removeMomentum in MDGenerator works properly or not.
 */
void InputFileTest::testRemoveMomentum() {

	std::vector<Component>& components = _domain->getComponents();

	ParticleContainer* particleContainer
		= initializeFromFile(ParticleContainerFactory::LinkedCell, "1clj-regular-2x2x3-removeMomentum.inp", 1.8);

	MDGenerator::removeMomentum(particleContainer, components);

	double mass;
	double mass_sum=0.;
	double momentum_sum[3] = {0., 0., 0.};

	Molecule* molecule = particleContainer->begin();
	while(molecule != particleContainer->end()){
		mass = components[molecule->componentid()].m();
		mass_sum = mass_sum + mass;
		momentum_sum[0] = momentum_sum[0] + mass * molecule->v(0);
		momentum_sum[1] = momentum_sum[1] + mass * molecule->v(1);
		momentum_sum[2] = momentum_sum[2] + mass * molecule->v(2);
		molecule = particleContainer->next();
	}

	ASSERT_DOUBLES_EQUAL(0., momentum_sum[0],1e-6);
	ASSERT_DOUBLES_EQUAL(0., momentum_sum[1],1e-6);
	ASSERT_DOUBLES_EQUAL(0., momentum_sum[2],1e-6);

	delete particleContainer;
}
#endif

