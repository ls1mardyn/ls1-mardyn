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
	TEST_METHOD(testExchangeMolecules);
	TEST_SUITE_END();

public:
	DomainDecompBaseTest();

	virtual ~DomainDecompBaseTest();

	void testExchangeMolecules();
};

#endif /* DOMAINDECOMPBASETEST_H_ */
