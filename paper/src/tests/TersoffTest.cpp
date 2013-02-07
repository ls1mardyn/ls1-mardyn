/*
 * TersoffTest.h
 *
 * @Date: 10.08.2012
 * @Author: Florian Groetzner
 */

#include "TersoffTest.h"
#include "molecules/Molecule.h"

TEST_SUITE_REGISTRATION(TersoffTest);

TersoffTest::TersoffTest() {

}

TersoffTest::~TersoffTest() {
}


void TersoffTest::testTersoffbij() {
	double test = 0.0;
	// Components for the molecules to test
	std::vector<Component> components;
	Component dummyComponent(0);

	// Add tersoff
	dummyComponent.addTersoff(0,0,0,0,51.212,12.742,0,0,1.95,0.15,38049,4.3484,-0.57058,0.72751,1.5724e-7);
	components.push_back(dummyComponent);

	Molecule a(0, 0, 1.0, 1.0, 1.0,0,0,0,0,0,0,0,0,0,0, &components);
	Molecule b(0, 0, 2.0, 2.0, 2.0,0,0,0,0,0,0,0,0,0,0, &components);

	// Initialize parameter list
	double params[15] = {0};
	a.tersoffParameters(params);

	double distance[3] = {1.0,0.0,0.0};

	test = Tersoffbij(&a, &b, params, distance, 1.0);

	a.addTersoffNeighbour(&b, false);
	b.addTersoffNeighbour(&a, false);


	//ASSERT_EQUAL(0ul, 1ul);
}
