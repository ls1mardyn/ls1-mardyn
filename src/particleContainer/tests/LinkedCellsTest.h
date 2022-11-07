/*
 * LinkedCellsTest.h
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#ifndef LINKEDCELLSTEST_H_
#define LINKEDCELLSTEST_H_

#include "particleContainer/tests/ParticleContainerTest.h"


#include "particleContainer/TraversalTuner.h"

class LinkedCellsTest: public ParticleContainerTest {

	TEST_SUITE(LinkedCellsTest);
	TEST_METHOD(testInsertion);
	TEST_METHOD(testMoleculeIteration);
	TEST_METHOD(testUpdateAndDeleteOuterParticles);
	TEST_METHOD(testUpdateAndDeleteOuterParticlesH2O);
	TEST_METHOD(testUpdateAndDeleteOuterParticles8Particles);
	TEST_METHOD(testMoleculeBeginNextEndDeleteCurrent);
	TEST_METHOD(testTraversalMethods);

	TEST_METHOD(testRegionIterator);
	TEST_METHOD(testRegionIteratorFile);

	TEST_METHOD(testCellBorderAndFlagManager);

#ifndef ENABLE_REDUCED_MEMORY_MODE
	TEST_METHOD(testFullShellMPIDirectPP);
	TEST_METHOD(testFullShellMPIDirect);

	TEST_METHOD(testHalfShellMPIDirectPP);
	TEST_METHOD(testHalfShellMPIDirect);
	TEST_METHOD(testHalfShellMPIIndirect);

	TEST_METHOD(testMidpointMPIDirectPP);
	TEST_METHOD(testMidpointMPIDirect);
	TEST_METHOD(testMidpointMPIIndirect);

	TEST_METHOD(testEighthShellMPIDirectPP);
#else
#pragma message "half and midpoint tests disabled for RMM"
#endif

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
	void testTraversalMethods();
	void testRegionIterator();
	void testRegionIteratorFile();
	void testGetHaloBoundaryParticlesDirection();

	void testHalfShell();

	void testFullShellMPIDirectPP();
	void testFullShellMPIDirect();

	void testHalfShellMPIDirectPP();
	void testHalfShellMPIDirect();
	void testHalfShellMPIIndirect();

	void testMidpoint();
	void testMidpointMPIDirectPP();
	void testMidpointMPIDirect();
	void testMidpointMPIIndirect();

	void testEighthShellMPIDirectPP();

	void testCellBorderAndFlagManager();

private:

	void doForceComparisonTest(std::string inputFile, TraversalTuner<ParticleCell>::traversalNames traversal, unsigned cellsInCutoff, std::string neighbourCommScheme, std::string commScheme);
};

#endif /* LINKEDCELLSTEST_H_ */
