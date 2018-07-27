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

	vector<int> v;
	int array[] = { 5, 4, 2, 3 };

	vector<int>::iterator it = v.begin();
	v.insert(it, array, array+4);
	it = v.begin();

	UnorderedVector::fastRemove(v, it);
	// 5 was erased, 3 was copied in its place, vector contains 3, 4, 2, iterator points at 3
	ASSERT_EQUAL(*it, 3);

	++it;
	UnorderedVector::fastRemove(v, it);
	// 4 was erased, 2 was copied in its place, vector contains 3, 2, iterator points at 2
	ASSERT_EQUAL(*it, 2);

	UnorderedVector::fastRemove(v, it);
	// 2 was erased, vector contains 3, iterator points at end()
	ASSERT_TRUE(it == v.end());

	it = v.begin();
	UnorderedVector::fastRemove(v, it);
	// 2 was erased, vector is empty and iterator points at end()
	ASSERT_TRUE(it == v.end());
	ASSERT_TRUE(v.empty());
}

void UnorderedVectorTest::testFastRemovalMoleculePointer() {
	using std::vector;

	vector<Molecule *> v;
	vector<Molecule*>::iterator it = v.begin();
	Component c(0);

	unsigned long ids[4] = { 5ul, 4ul, 2ul, 3ul };

	for(int i = 0; i < 4; ++i) {
		v.push_back(
				new Molecule(ids[i], &c, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
						0., 0., 0.));
	}

	// erase 5:

	it = v.begin();
	ASSERT_EQUAL((*it)->id(), 5ul);

	delete *it;
	UnorderedVector::fastRemove(v, it);
	// 5 was erased, 3 was copied in its place, vector contains 3, 4, 2, iterator points at 3
	ASSERT_EQUAL((*it)->id(), 3ul);


	// erase 4:

	++it;
	ASSERT_EQUAL((*it)->id(), 4ul);
	delete *it;
	UnorderedVector::fastRemove(v, it);
	// 4 was erased, 2 was copied in its place, vector contains 3, 2, iterator points at 2
	ASSERT_EQUAL((*it)->id(), 2ul);


	// erase 2:

	delete *it;
	UnorderedVector::fastRemove(v, it);
	// 2 was erased, vector contains 3, iterator points at end()
	ASSERT_TRUE(it == v.end());


	// erase 3:

	it = v.begin();
	delete *it;
	UnorderedVector::fastRemove(v, it);
	// 3 was erased, vector is empty and iterator points at end()
	ASSERT_TRUE(it == v.end());
	ASSERT_TRUE(v.empty());

}
