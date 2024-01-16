/*
 * UnorderedVectorTest.cpp
 *
 *  Created on: Nov 19, 2015
 *      Author: tchipevn
 */

#include "UnorderedVectorTest.h"
#include "../UnorderedVector.h"
#include "molecules/Molecule.h"
#include <vector>

TEST_SUITE_REGISTRATION(UnorderedVectorTest);

void UnorderedVectorTest::testFastRemovalInt() {
	using std::vector;

	std::vector<int> v;
	int array[] = { 5, 4, 2, 3 };

	v.insert(v.begin(), array, array+4);
	unsigned int index = 0;

	UnorderedVector::fastRemove(v, index);
	// 5 was erased, 3 was copied in its place, vector contains 3, 4, 2, iterator points at 3
	ASSERT_EQUAL(v[index], 3);

	++index;
	UnorderedVector::fastRemove(v, index);
	// 4 was erased, 2 was copied in its place, vector contains 3, 2, iterator points at 2
	ASSERT_EQUAL(v[index], 2);

	UnorderedVector::fastRemove(v, index);
	// 2 was erased, vector contains 3, iterator points at end()
	ASSERT_TRUE(index == v.size());

	index = 0;
	UnorderedVector::fastRemove(v, index);
	// 2 was erased, vector is empty and iterator points at end()
	ASSERT_TRUE(index == v.size());
	ASSERT_TRUE(v.empty());
}

void UnorderedVectorTest::testFastRemovalMoleculePointer() {
	using std::vector;

	std::vector<Molecule *> v;
	unsigned int index;
	Component c(0);

	unsigned long ids[4] = { 5ul, 4ul, 2ul, 3ul };

	for(int i = 0; i < 4; ++i) {
		v.push_back(
				new Molecule(ids[i], &c, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
						0., 0., 0.));
	}

	// erase 5:

	index = 0;
	ASSERT_EQUAL(v[index]->getID(), 5ul);

	delete v[index];
	UnorderedVector::fastRemove(v, index);
	// 5 was erased, 3 was copied in its place, vector contains 3, 4, 2, iterator points at 3
	ASSERT_EQUAL((v[index])->getID(), 3ul);


	// erase 4:

	++index;
	ASSERT_EQUAL((v[index])->getID(), 4ul);
	delete v[index];
	UnorderedVector::fastRemove(v, index);
	// 4 was erased, 2 was copied in its place, vector contains 3, 2, iterator points at 2
	ASSERT_EQUAL((v[index])->getID(), 2ul);


	// erase 2:

	delete v[index];
	UnorderedVector::fastRemove(v, index);
	// 2 was erased, vector contains 3, iterator points at end()
	ASSERT_TRUE(index == v.size());


	// erase 3:

	index = 0;
	delete v[index];
	UnorderedVector::fastRemove(v, index);
	// 3 was erased, vector is empty and iterator points at end()
	ASSERT_TRUE(index == v.size());
	ASSERT_TRUE(v.empty());

}
