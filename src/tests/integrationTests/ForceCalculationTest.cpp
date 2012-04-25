/*
 * ForceCalculationTest.cpp
 *
 * @Date: 22.04.2012
 * @Author: eckhardw
 */

#include "ForceCalculationTest.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

TEST_SUITE_REGISTRATION(ForceCalculationTest);

ForceCalculationTest::ForceCalculationTest() {
}

ForceCalculationTest::~ForceCalculationTest() {
}

void ForceCalculationTest::testForcePotentialCalculation() {
	// U (r_ij) = 4 epsilon * ( (sigma / r_ij)^12 - (sigma / r_ij)^6 )
	// F (r_ij) = 24 epsilon * 1/r ( (sigma / r_ij)^6 - 2 * (sigma / r_ij)^12 )

	double forces[4][3] = { { -24, -24, 0 },
	                        {  24, -24, 0 },
	                        { -24,  24, 0 },
	                        {  24,  24, 0 }};

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "1clj-regular-2x2.inp", 1.1);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-8);

	ParticlePairs2PotForceAdapter forceAdapter(*_domain);
	container->traversePairs(&forceAdapter);

	for (Molecule* m = container->begin(); m != container->end(); m = container->next()) {
		m->calcFM();
	}

	for (Molecule* m = container->begin(); m != container->end(); m = container->next()) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m->id() << " index i="<< i << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), forces[m->id()-1][i], m->F(i), 1e-8);
		}
	}

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(96, _domain->getLocalVirial(), 1e-8);

	delete container;
}
