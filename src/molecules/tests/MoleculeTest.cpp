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
    decltype(std::declval<Component>().ID()) const id = 0;
    Component dummyComponent(0);
	Molecule a(42,                //id
			&dummyComponent,      //component id
			1.0,2.0,3.0,          //r
			4.0,5.0,6.0,          //v
			7.0,8.0,9.0,10.0,     //q
			11.0,12.0,13.0);      //L
    std::vector<char> buffer(a.serializedSize());
	auto lastElement = a.serialize(buffer.begin());

    //these unit test will all break if the type of the molecule interface data changes
	ASSERT_TRUE(!(lastElement-buffer.end()));
    ASSERT_TRUE(*reinterpret_cast<unsigned long*>(&buffer[0]) == 42);
    ASSERT_TRUE(*reinterpret_cast<decltype(id)*>(&buffer[8]) == id+1);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[12]) == 1.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[20]) == 2.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[28]) == 3.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[36]) == 4.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[44]) == 5.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[52]) == 6.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[60]) == 7.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[68]) == 8.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[76]) == 9.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[84]) == 10.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[92]) == 11.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[100]) == 12.0);
    ASSERT_TRUE(*reinterpret_cast<double*>(&buffer[108]) == 13.0);
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
