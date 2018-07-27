/*
 * ParticleCellRMMTest.h
 *
 *  Created on: 12 Oct 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_TESTS_PARTICLECELLRMMTEST_H_
#define SRC_PARTICLECONTAINER_TESTS_PARTICLECELLRMMTEST_H_

#include "utils/Testing.h"

class ParticleCellRMMTest : public utils::Test {

	TEST_SUITE(ParticleCellRMMTest);
	TEST_METHOD(testSizeOfIs64);
	TEST_SUITE_END();

public:
	ParticleCellRMMTest();
	virtual ~ParticleCellRMMTest();

	void testSizeOfIs64();
};

#endif /* SRC_PARTICLECONTAINER_TESTS_PARTICLECELLRMMTEST_H_ */
