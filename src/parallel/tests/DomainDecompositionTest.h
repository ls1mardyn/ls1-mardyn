/*
 * DomainDecompositionTest.h
 *
 * @Date: 10.05.2012
 * @Author: eckhardw
 */

#ifndef DOMAINDECOMPOSITIONTEST_H_
#define DOMAINDECOMPOSITIONTEST_H_

#include "utils/TestWithSimulationSetup.h"

class DomainDecompositionTest: public utils::TestWithSimulationSetup {

	TEST_SUITE(DomainDecompositionTest);
	TEST_METHOD(testNoDuplicatedParticles);
	TEST_METHOD(testNoLostParticles);
	TEST_METHOD(testExchangeMolecules1Proc);
	TEST_SUITE_END();

public:

	DomainDecompositionTest();

	~DomainDecompositionTest() override;

	void testNoDuplicatedParticles();
	void testNoLostParticles();
	/**
	 * Test the particle exchange if running with 1 process.
	 */
	void testExchangeMolecules1Proc();
private:
	void testNoDuplicatedParticlesFilename(const char * filename, double cutoff);
	void testNoLostParticlesFilename(const char * filename, double cutoff);
};

#endif /* DOMAINDECOMPOSITIONTEST_H_ */
