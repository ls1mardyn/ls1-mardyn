/*
 * CompareFMMContainersTest.cpp
 *
 *  Created on: Nov 19, 2015
 *      Author: tchipevn
 */

#include "CompareFMMContainersTest.h"
#include "molecules/Molecule.h"
#include "bhfmm/FastMultipoleMethod.h"

TEST_SUITE_REGISTRATION(CompareFMMContainersTest);

CompareFMMContainersTest::CompareFMMContainersTest() {
	// TODO Auto-generated constructor stub

}

CompareFMMContainersTest::~CompareFMMContainersTest() {
	// TODO Auto-generated destructor stub
}

void CompareFMMContainersTest::compare(double cutoffRadius, bool periodic) {

	double globalDomainLength[3] = {8., 8., 8.};
	double bBoxMin[3] = {0., 0., 0.};
	double bBoxMax[3] = {8., 8., 8.};
	double LJCellLength[3] = {cutoffRadius, cutoffRadius, cutoffRadius};
	unsigned LJSubdivisionFactor = 1;
	int orderOfExpansions = 2;

	// Uniform container
	ParticleContainer * LCUniform = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", cutoffRadius);
	bool adaptiveArg = false;

	bhfmm::FastMultipoleMethod uniform;
	uniform.setParameters(LJSubdivisionFactor, orderOfExpansions, periodic, adaptiveArg);
	uniform.init(globalDomainLength, bBoxMin, bBoxMax, LJCellLength);

	uniform.computeElectrostatics(LCUniform);
	for (Molecule* m = LCUniform->begin(); m != LCUniform->end(); m = LCUniform->next()) {
		m->calcFM();
	}

	// reset variables, which are not visible here
	tearDown();
	setUp();

	// Adaptive container
	ParticleContainer * LCAdaptive = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", cutoffRadius);
	adaptiveArg = true;

	bhfmm::FastMultipoleMethod adaptive;
	adaptive.setParameters(LJSubdivisionFactor, orderOfExpansions, periodic, adaptiveArg);
	adaptive.init(globalDomainLength, bBoxMin, bBoxMax, LJCellLength);

	adaptive.computeElectrostatics(LCAdaptive);
	for (Molecule* m = LCAdaptive->begin(); m != LCAdaptive->end(); m = LCAdaptive->next()) {
		m->calcFM();
	}

	// traverse molecules and compare forces
	Molecule * itUniform = LCUniform->begin();
	Molecule * itAdaptive = LCAdaptive->begin();
	for(; itUniform != LCUniform->end() and itAdaptive != LCAdaptive->end(); itUniform = LCUniform->next(), itAdaptive = LCAdaptive->next()) {
		ASSERT_DOUBLES_EQUAL_MSG("Force component x should be equal", itUniform->F(0), itAdaptive->F(0), 1e-12);
		ASSERT_DOUBLES_EQUAL_MSG("Force component y should be equal", itUniform->F(1), itAdaptive->F(1), 1e-12);
		ASSERT_DOUBLES_EQUAL_MSG("Force component z should be equal", itUniform->F(2), itAdaptive->F(2), 1e-12);
	}

	delete LCUniform;
	delete LCAdaptive;
}

void CompareFMMContainersTest::compareAtRadius1() {
	compare(1.0);
}

void CompareFMMContainersTest::compareAtRadius2() {
	compare(2.0);
}

void CompareFMMContainersTest::compareAtRadius4() {
	compare(4.0);
}

void CompareFMMContainersTest::compareAtRadius1WithoutPeriodicBC() {
	compare(1.0, false);
}

void CompareFMMContainersTest::compareAtRadius2WithoutPeriodicBC() {
	compare(2.0, false);
}

void CompareFMMContainersTest::compareAtRadius4WithoutPeriodicBC() {
	compare(4.0, false);
}
