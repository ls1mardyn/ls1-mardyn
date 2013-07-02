/*
 * MoleculeTest.cpp
 *
 * @Date: 04.02.2011
 * @Author: eckhardw
 */

#include "MoleculeTest.h"
#include "molecules/Molecule.h"

TEST_SUITE_REGISTRATION(MoleculeTest);

MoleculeTest::MoleculeTest() {

}

MoleculeTest::~MoleculeTest() {
}


void MoleculeTest::testIsLessThan() {
	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,0,0,0,0,false);
	components.push_back(dummyComponent);

	Molecule a(0, 0, 1.0, 1.0, 1.0,0,0,0,0,0,0,0,0,0,0, &components);
	Molecule b(0, 0, 2.0, 2.0, 2.0,0,0,0,0,0,0,0,0,0,0, &components);

	ASSERT_TRUE(a.isLessThan(b));
	ASSERT_TRUE(!b.isLessThan(a));

	a.setr(2, 3.0);

	ASSERT_TRUE(!a.isLessThan(b));
	ASSERT_TRUE(b.isLessThan(a));

	a.setr(2, 2.0);
	ASSERT_TRUE(a.isLessThan(b));
	ASSERT_TRUE(!b.isLessThan(a));

	a.setr(1, 3.0);

	ASSERT_TRUE(!a.isLessThan(b));
	ASSERT_TRUE(b.isLessThan(a));

	a.setr(1, 2.0);
	ASSERT_TRUE(a.isLessThan(b));
	ASSERT_TRUE(!b.isLessThan(a));

	a.setr(0, 3.0);

	ASSERT_TRUE(!a.isLessThan(b));
	ASSERT_TRUE(b.isLessThan(a));
}
