/*
 * AdaptiveSubCellsTest.h
 *
 * @Date: 09.05.2012
 * @Author: eckhardw
 */

#ifndef ADAPTIVESUBCELLSTEST_H_
#define ADAPTIVESUBCELLSTEST_H_

#include "particleContainer/tests/ParticleContainerTest.h"
#include "particleContainer/AdaptiveSubCells.h"

class AdaptiveSubCellsTest : public ParticleContainerTest{

	TEST_SUITE(AdaptiveSubCellsTest);
	TEST_METHOD(testInsertion);
	TEST_METHOD(testMoleculeIteration);
	TEST_METHOD(testUpdateAndDeleteOuterParticles);
	TEST_SUITE_END();

public:

	AdaptiveSubCellsTest();

	virtual ~AdaptiveSubCellsTest();


	void testInsertion() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		AdaptiveSubCells container(boundings_min, boundings_max, 2.5, 2.5);
		this->ParticleContainerTest::testInsertion(&container);
	}

	void testMoleculeIteration() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		AdaptiveSubCells container(boundings_min, boundings_max, 2.5, 2.5);
		this->ParticleContainerTest::testMoleculeIteration(&container);
	}

	void testUpdateAndDeleteOuterParticles() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		AdaptiveSubCells container(boundings_min, boundings_max, 2.5, 2.5);
		this->ParticleContainerTest::testUpdateAndDeleteOuterParticles(&container);
	}
};

#endif /* ADAPTIVESUBCELLSTEST_H_ */
