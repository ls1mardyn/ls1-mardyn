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
		
	int numTests = 5;
	// Test results
	//double expected[numTests] = {1.0};
	double test = 0.0;

	// Components for the molecules to test
	std::vector<Component> components;
	Component dummyComponent(0);

	// Add tersoff
	dummyComponent.addTersoff(0,0,0,0,51.212,12.742,0,0,1.95,0.15,38049,4.3484,-0.57058,0.72751,1.5724e-7);
	components.push_back(dummyComponent);

	Molecule a(0, 0, 1.0, 1.0, 1.0,0,0,0,0,0,0,0,0,0,0, &components);

	// Initialize parameter list
	double params[15] = {0};
	a.tersoffParameters(params);
	
	// Distance vector
	double distance[3] = {1.0,0.0,0.0};

	double dis;	
	for (int i = 0; i < numTests; i++){	
		dis = (0.1*i) - 10.0;

		// Distance vector
		double distance[3] = {dis,0.0,0.0};

		Molecule b(1, 1, 0.0, 0.0, dis,0,0,0,0,0,0,0,0,0,0, &components);
		
		// Add neighbours
		a.addTersoffNeighbour(&b, true);
		b.addTersoffNeighbour(&a, true);

		test = Tersoffbij(&a, &b, params, distance, dis);
	
		//ASSERT_DOUBLES_EQUAL(*expected[i], test, 1e-6);
	}

}


void TersoffTest::tersoffUIJplusUJI() {
	
	int numTests = 5;
	// Test results
	//double expected[numTests] = {1.0};
	double test = 0.0;

	// Components for the molecules to test
	std::vector<Component> components;
	Component dummyComponent(0);

	// Add tersoff
	dummyComponent.addTersoff(0,0,0,0,51.212,12.742,0,0,1.95,0.15,38049,4.3484,-0.57058,0.72751,1.5724e-7);
	components.push_back(dummyComponent);

	Molecule a(0, 0, 1.0, 1.0, 1.0,0,0,0,0,0,0,0,0,0,0, &components);

	// Initialize parameter list
	double params[15] = {0};
	a.tersoffParameters(params);

	double upotTersoff = 1.0;

	double dis;	
	for (int i = 0; i < numTests; i++){
		dis = (0.1*i) - 10.0;

		Molecule b(1, 1, 0.0, 0.0, dis,0,0,0,0,0,0,0,0,0,0, &components);
		
		// Add neighbours
		a.addTersoffNeighbour(&b, true);
		b.addTersoffNeighbour(&a, true);

		// Distance vector
		double distance[3] = {dis,0.0,0.0};

		test = TersoffUIJplusUJI(&a, &b, params, distance, dis, upotTersoff);
	
		//ASSERT_DOUBLES_EQUAL(*expected[i], test, 1e-6);
	}

}

void TersoffTest::tersoffUIJattr() {
	
	int numTests = 5;
	// Test results
	//double expected[numTests] = {1.0};
	double test = 0.0;

	// Components for the molecules to test
	std::vector<Component> components;
	Component dummyComponent(0);

	// Add tersoff
	dummyComponent.addTersoff(0,0,0,0,51.212,12.742,0,0,1.95,0.15,38049,4.3484,-0.57058,0.72751,1.5724e-7);
	components.push_back(dummyComponent);

	Molecule a(0, 0, 1.0, 1.0, 1.0,0,0,0,0,0,0,0,0,0,0, &components);

	// Initialize parameter list
	double params[15] = {0};
	a.tersoffParameters(params);
	
	double dis;	
	for (int i = 0; i < numTests; i++){
		dis = (0.1*i) - 10.0;

		Molecule b(1, 1, 0.0, 0.0, dis,0,0,0,0,0,0,0,0,0,0, &components);
		
		// Add neighbours
		a.addTersoffNeighbour(&b, true);
		b.addTersoffNeighbour(&a, true);

		// Distance vector
		double distance[3] = {dis,0.0,0.0};

		test = TersoffUIJattr(&a, &b, params, distance, dis);
	
		//ASSERT_DOUBLES_EQUAL(*expected[i], test, 1e-6);
	}

}

void TersoffTest::tersoffPotential() {
	
	int numTests = 5;
	// Test results
	double expected[5] = {1.7674e31, 6.9530e29, 2.7354e28, 1.0761e27, 4.2335e25};
	double test = 0.0;

	// Components for the molecules to test
	std::vector<Component> components;
	Component dummyComponent(0);

	// Add tersoff
	dummyComponent.addTersoff(0,0,0,0,51.212,12.742,0,0,1.95,0.15,38049,4.3484,-0.57058,0.72751,1.5724e-7);
	components.push_back(dummyComponent);

	Molecule a(0, 0, 0.0, 0.0, 0.0,0,0,0,0,0,0,0,0,0,0, &components);
	
	// Initialize parameter list
	double params[15] = {0};
	a.tersoffParameters(params);
	
	double upotTersoff = 1.0;
	
	for (int i = 0; i < numTests; i++){
		Molecule b(1, 1, 0.0, 0.0, (0.1*i) - 10.0,0,0,0,0,0,0,0,0,0,0, &components);
		
		// Add neighbours
		a.addTersoffNeighbour(&b, true);
		b.addTersoffNeighbour(&a, true);

		test = TersoffPotential(&a, params, upotTersoff);
		
		double result = expected[i];
		// ASSERT_DOUBLES_EQUAL(result, test, 1e-1);
	}

}



