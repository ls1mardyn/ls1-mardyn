/*
 * DttNodeTest.h
 *
 *  Created on: Nov 19, 2015
 *      Author: tchipevn
 */

#ifndef SRC_BHFMM_CONTAINERS_TESTS_DTTNODETEST_H_
#define SRC_BHFMM_CONTAINERS_TESTS_DTTNODETEST_H_

#include "utils/TestWithSimulationSetup.h"
#include "bhfmm/containers/DttNode.h"

class DttNodeTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(DttNodeTest);

	TEST_METHOD(testDepthAtRadius1);
	TEST_METHOD(testDepthAtRadius2);
	TEST_METHOD(testDepthAtRadius4);

	TEST_METHOD(testDivideParticles);

	TEST_METHOD(testSoAConvertions);

	TEST_METHOD(testUpwardDownwardWithNoInteraction);

	TEST_SUITE_END();

public:
	DttNodeTest();
	virtual ~DttNodeTest();

	void testDepthAtRadius1();
	void testDepthAtRadius2();
	void testDepthAtRadius4();

	void testDivideParticles();

	void testSoAConvertions();
	void testUpwardDownwardWithNoInteraction();
private:
	void testDepth(double cutoffRadius);
};

#endif /* SRC_BHFMM_CONTAINERS_TESTS_DTTNODETEST_H_ */
