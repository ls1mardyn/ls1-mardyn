/*
 * DomainDecompBaseTest.h
 *
 * @Date: 09.05.2012
 * @Author: eckhardw
 */

#ifndef DOMAINDECOMPBASETEST_H_
#define DOMAINDECOMPBASETEST_H_

#include "utils/TestWithSimulationSetup.h"

class DomainDecompBaseTest: public utils::TestWithSimulationSetup {

	TEST_SUITE(DomainDecompBaseTest);
	TEST_METHOD(testNoDuplicatedParticles);
	TEST_METHOD(testNoLostParticles);
	TEST_METHOD(testExchangeMoleculesSimple);
	TEST_METHOD(testExchangeMolecules);
	TEST_SUITE_END();

public:
	DomainDecompBaseTest();

	virtual ~DomainDecompBaseTest();

	void testNoDuplicatedParticles();

	void testNoLostParticles();
	void testExchangeMoleculesSimple();
	void testExchangeMolecules();

private:
	void testNoDuplicatedParticlesFilename(const char * filename, double cutoff);
	void testNoLostParticlesFilename(const char * filename, double cutoff);
};

#endif /* DOMAINDECOMPBASETEST_H_ */
