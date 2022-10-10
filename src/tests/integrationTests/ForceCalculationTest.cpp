/*
 * ForceCalculationTest.cpp
 *
 * @Date: 22.04.2012
 * @Author: eckhardw
 */

#include "ForceCalculationTest.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"

TEST_SUITE_REGISTRATION(ForceCalculationTest);

ForceCalculationTest::ForceCalculationTest() {
}

ForceCalculationTest::~ForceCalculationTest() {
}

void ForceCalculationTest::testForcePotentialCalculationU0() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "ForceCalculationTest::testForcePotentialCalculationU0: SKIPPED (required exactly 1 process but was run with " <<  _domainDecomposition->getNumProcs() << " processes)" << std::endl;
		return;
	}

	// U (r_ij) = 4 epsilon * ( (sigma / r_ij)^12 - (sigma / r_ij)^6 )
	// F (r_ij) = 24 epsilon * 1/r ( (sigma / r_ij)^6 - 2 * (sigma / r_ij)^12 )

	double forces[4][3] = { { -24, -24, 0 },
	                        {  24, -24, 0 },
	                        { -24,  24, 0 },
	                        {  24,  24, 0 }};

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "ForceCalculationTestU0.inp", 1.1);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-8);

	ParticlePairs2PotForceAdapter forceAdapter(*_domain);
	LegacyCellProcessor cellProcessor( 1.1, 1.1, &forceAdapter);
	container->traverseCells(cellProcessor);

	for (auto m = container->iterator(ParticleIterator::ALL_CELLS); m.isValid(); ++m) {
		m->calcFM();
	}

	for (auto m = container->iterator(ParticleIterator::ALL_CELLS); m.isValid(); ++m) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m->getID() << " index i="<< i << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), forces[m->getID()-1][i], m->F(i), 1e-8);
		}
	}

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(96, _domain->getLocalVirial(), 1e-8);

	delete container;
}

void ForceCalculationTest::testForcePotentialCalculationF0() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "ForceCalculationTest::testForcePotentialCalculationF0: SKIPPED (required exactly 1 process but was run with " <<  _domainDecomposition->getNumProcs() << " processes)" << std::endl;
		return;
	}
#if defined(MARDYN_DPDP)
	double tolerance_force = 1e-7;
	double tolerance_virial = 1e-6;
#else
	double tolerance_force = 1e-6;
	double tolerance_virial = 1e-5;
#endif


	// U (r_ij) = 4 epsilon * ( (sigma / r_ij)^12 - (sigma / r_ij)^6 )
	// F (r_ij) = 24 epsilon * 1/r ( (sigma / r_ij)^6 - 2 * (sigma / r_ij)^12 )
	// -> potential U = -1 per particle pair.
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "ForceCalculationTestF0.inp", 1.3);

    ASSERT_EQUAL(4ul, container->getNumberOfParticles());
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-8);

	ParticlePairs2PotForceAdapter forceAdapter(*_domain);
	LegacyCellProcessor cellProcessor( 1.3, 1.3, &forceAdapter);
	container->traverseCells(cellProcessor);

	for (auto m = container->iterator(ParticleIterator::ALL_CELLS); m.isValid(); ++m) {
		m->calcFM();
	}

	for (auto m = container->iterator(ParticleIterator::ALL_CELLS); m.isValid(); ++m) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m->getID() << " index i="<< i << " F[i]=" << m->F(i) << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), 0.0, m->F(i), tolerance_force);
		}
	}

	ASSERT_DOUBLES_EQUAL(-4, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), tolerance_virial);

	delete container;
}
