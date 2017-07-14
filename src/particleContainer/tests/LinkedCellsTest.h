/*
 * LinkedCellsTest.h
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#ifndef LINKEDCELLSTEST_H_
#define LINKEDCELLSTEST_H_

#include "particleContainer/tests/ParticleContainerTest.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/TraversalTuner.h"

class LinkedCellsTest: public ParticleContainerTest {

	TEST_SUITE(LinkedCellsTest);
	TEST_METHOD(testInsertion);
	TEST_METHOD(testMoleculeIteration);
	TEST_METHOD(testUpdateAndDeleteOuterParticles);
	TEST_METHOD(testUpdateAndDeleteOuterParticlesH2O);
	TEST_METHOD(testUpdateAndDeleteOuterParticles8Particles);
	TEST_METHOD(testMoleculeBeginNextEndDeleteCurrent);
	TEST_METHOD(testParticleIteratorBeginNextEndParticleIteratorSequential);
	//TEST_METHOD(testHalfShell);
	TEST_METHOD(testHalfShellMPIDirect);
	TEST_METHOD(testHalfShellMPIIndirect);
//	TEST_METHOD(testGetHaloBoundaryParticlesDirection);
	TEST_METHOD(testTraversalMethods);
	//TEST_METHOD(testMidpoint);
	TEST_METHOD(testMidpointMPI);
	TEST_SUITE_END();

public:

	LinkedCellsTest();

	virtual ~LinkedCellsTest();

	void testInsertion() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		LinkedCells container(boundings_min, boundings_max, 2.5);
		this->ParticleContainerTest::testInsertion(&container);
	}

	void testMoleculeIteration() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		LinkedCells container(boundings_min, boundings_max, 2.5);
		this->ParticleContainerTest::testMoleculeIteration(&container);
	}

	void testUpdateAndDeleteOuterParticles() {
		double boundings_min[] = {0, 0, 0};
		double boundings_max[] = {10.0, 10.0, 10.0 };
		LinkedCells container(boundings_min, boundings_max, 2.5);
		this->ParticleContainerTest::testUpdateAndDeleteOuterParticles(&container);
	}

	void testUpdateAndDeleteOuterParticlesFilename(const char * filename, double cutoff);
	void testUpdateAndDeleteOuterParticlesH2O();
	void testUpdateAndDeleteOuterParticles8Particles();
	void testMoleculeBeginNextEndDeleteCurrent();
	void testParticleIteratorBeginNextEndParticleIteratorSequential();
	void testTraversalMethods();
//	void testGetHaloBoundaryParticlesDirection();

	//void testHalfShell();
	void testHalfShellMPIDirect();
	void testHalfShellMPIIndirect();

	//void testMidpoint();
	void testMidpointMPI();

private:

	void doForceComparisonTest(std::string inputFile, TraversalTuner::traversalNames traversal, std::string neighbourCommScheme, std::string commScheme);
};

#endif /* LINKEDCELLSTEST_H_ */
