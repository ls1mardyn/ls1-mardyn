/*
 * CompareFMMContainersTest.h
 *
 *  Created on: Nov 19, 2015
 *      Author: tchipevn
 */

#ifndef SRC_BHFMM_TESTS_COMPAREFMMCONTAINERSTEST_H_
#define SRC_BHFMM_TESTS_COMPAREFMMCONTAINERSTEST_H_

#include "utils/TestWithSimulationSetup.h"

class CompareFMMContainersTest: public utils::TestWithSimulationSetup {
	TEST_SUITE(CompareFMMContainersTest);

	TEST_METHOD(compareAtRadius4WithoutPeriodicBC);
	TEST_METHOD(compareAtRadius2WithoutPeriodicBC);
	TEST_METHOD(compareAtRadius1WithoutPeriodicBC);

	/*TEST_METHOD(compareAtRadius4);
	TEST_METHOD(compareAtRadius2);
	TEST_METHOD(compareAtRadius1);*/

	TEST_SUITE_END();

public:
	CompareFMMContainersTest();
	virtual ~CompareFMMContainersTest();

	void compareAtRadius1();
	void compareAtRadius2();
	void compareAtRadius4();
	void compareAtRadius1WithoutPeriodicBC();
	void compareAtRadius2WithoutPeriodicBC();
	void compareAtRadius4WithoutPeriodicBC();

private:
	void compare(double cutoffRadius, bool periodic = true);

};

#endif /* SRC_BHFMM_TESTS_COMPAREFMMCONTAINERSTEST_H_ */
