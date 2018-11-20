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

void MoleculeTest::testSerialize() {
	Component dummyComponent(0);
	Molecule a(0,                 //id
			&dummyComponent,      //component id (address)???
			1.0,1.0,1.0,          //r
			0.0,0.0,0.0,          //v
			0.0,0.0,0.0,0.0,      //q
			0,0,0);               //L
	
	std::vector<char> buffer(a.serializedSize());
	auto lastElement = a.serialize(buffer.begin());
	// ASSERT_TRUE(!(lastElement-buffer.end()));
}

void MoleculeTest::testIsLessThan() {
	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,0,0,0,0,false);
	components.push_back(dummyComponent);

	Molecule a(0, &components[0], 1.0, 1.0, 1.0,0,0,0,0,0,0,0,0,0,0);
	Molecule b(0, &components[0], 2.0, 2.0, 2.0,0,0,0,0,0,0,0,0,0,0);

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
