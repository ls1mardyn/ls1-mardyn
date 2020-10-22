/*
 * KDNodeTest.cpp
 *
 *  Created on: Feb 24, 2012
 *      Author: eckhardw
 */

#include "KDNodeTest.h"
#include "parallel/KDNode.h"
#include <string>

TEST_SUITE_REGISTRATION(KDNodeTest);

KDNodeTest::KDNodeTest() {
}

KDNodeTest::~KDNodeTest() {
}


void KDNodeTest::testEqual() {
	int lowerEnd[] = {0, 0, 0};
	int upperEnd[] = {3, 3, 3};
	bool coversAll[] = {true, true, true};

	KDNode node1(5, lowerEnd, upperEnd, 41, 0, coversAll, 0);
	KDNode node2(5, lowerEnd, upperEnd, 41, 0, coversAll, 0);
	ASSERT_TRUE(node1.equals(node2));

	KDNode node3(7, lowerEnd, upperEnd, 41, 0, coversAll, 0);
	ASSERT_TRUE(! node1.equals(node3));

	KDNode* node4 = new KDNode(7, lowerEnd, upperEnd, 41, 0, coversAll, 0);
	node1._child1 = node4;
	KDNode* node5 = new KDNode(7, lowerEnd, upperEnd, 41, 0, coversAll, 0);
	node2._child1 = node5;
	// node4 and node5 are deleted by their root nodes...
	ASSERT_TRUE(node1.equals(node2));
}

void KDNodeTest::testSplit() {
	int lowerEnd[] = {0, 0, 0};
	int upperEnd[] = {7, 3, 3};
	bool coversAll[] = {true, true, true};

	KDNode root(7, lowerEnd, upperEnd, 1, 0, coversAll, 0);
	root.split(1, 1, 3);

	KDNode resultRoot(7, lowerEnd, upperEnd, 1, 0, coversAll, 0);
	int lower1[] = {0, 0, 0};
	int upper1[] = {7, 1, 3};
	bool childCoversAll[] = {true, false, true};
	KDNode* resultChild1 = new KDNode(3, lower1, upper1, 2, 0, childCoversAll, 1);

	int lower2[] = {0, 2, 0};
	int upper2[] = {7, 3, 3};
	KDNode* resultChild2 = new KDNode(4, lower2, upper2, 7, 3, childCoversAll, 1);
	resultRoot._child1 = resultChild1;
	resultRoot._child2 = resultChild2;

	ASSERT_TRUE(root.equals(resultRoot));
}

void KDNodeTest::testBuildKDTree() {

	int lowerEnd[] = {0, 0, 0};
	int upperEnd[] = {7, 3, 3};
	bool coversAll[] = {true, true, true};

	KDNode node1(1, lowerEnd, upperEnd, 1, 0, coversAll, 0);
	node1.buildKDTree();
	KDNode result(1, lowerEnd, upperEnd, 1, 0, coversAll, 0);
	ASSERT_TRUE(node1.equals(result));


	KDNode root(2, lowerEnd, upperEnd, 0, 0, coversAll, 0);
	root.buildKDTree();

	KDNode resultRoot(2, lowerEnd, upperEnd, 0, 0, coversAll, 0);
	int lower1[] = {0, 0, 0};
	int upper1[] = {3, 3, 3};
	bool childCoversAll[] = {false, true, true};
	KDNode* resultChild1 = new KDNode(1, lower1, upper1, 1, 0, childCoversAll, 1);

	int lower2[] = {4, 0, 0};
	int upper2[] = {7, 3, 3};
	KDNode* resultChild2 = new KDNode(1, lower2, upper2, 2, 1, childCoversAll, 1);
	resultRoot._child1 = resultChild1;
	resultRoot._child2 = resultChild2;

	ASSERT_TRUE(root.equals(resultRoot));
}

void KDNodeTest::testFindAreaForProcess() {
	int lowerEnd[] = {0, 0, 0};
	int upperEnd[] = {7, 3, 3};
	bool coversAll[] = {true, true, true};
	KDNode resultRoot(3, lowerEnd, upperEnd, 1, 0, coversAll, 0);

	int lower1[] = {0, 0, 0};
	int upper1[] = {3, 3, 3};
	bool childCoversAll[] = {false, true, true};
	KDNode* resultChild1 = new KDNode(2, lower1, upper1, 1, 0, childCoversAll, 1);
	int lower2[] = {4, 0, 0};
	int upper2[] = {7, 3, 3};
	KDNode* resultChild2 = new KDNode(1, lower2, upper2, 1, 2, childCoversAll, 1);

	resultRoot._child1 = resultChild1;
	resultRoot._child2 = resultChild2;

	int lower3[] = {0, 0, 0};
	int upper3[] = {1, 3, 3};
	KDNode* resultChild3 = new KDNode(1, lower3, upper3, 1, 0, childCoversAll, 2);

	int lower4[] = {2, 0, 0};
	int upper4[] = {3, 3, 3};
	KDNode* resultChild4 = new KDNode(1, lower4, upper4, 1, 1, childCoversAll, 2);

	resultRoot._child1->_child1 = resultChild3;
	resultRoot._child1->_child2 = resultChild4;

	ASSERT_EQUAL(resultRoot.findAreaForProcess(0), resultChild3);
	ASSERT_EQUAL(resultRoot.findAreaForProcess(1), resultChild4);
	ASSERT_EQUAL(resultRoot.findAreaForProcess(2), resultChild2);
	ASSERT_EQUAL(resultRoot.findAreaForProcess(41), (KDNode*) NULL);
}

void KDNodeTest::testGetMPIKDNode() {
	const int lowerCorner[3] = {0,0,0};
	const int higherCorner[3] = {3,3,3};
	std::bitset<3> coversWholeDomain;
	for (int i = 0; i < 3; i++) {
		coversWholeDomain[i] = true;
	}

	MPIKDNodePacked mpiNode(coversWholeDomain, 1, lowerCorner, higherCorner,
			7, 0, -1, -1, 1, 0.0, 0.0, 2.0, 3.0, 2);

	for (int i = 0; i < 3; i++) {
		ASSERT_EQUAL(mpiNode.getCoversWholeDomain(i), true);
		ASSERT_EQUAL(mpiNode.getLowCorner(i), 0);
		ASSERT_EQUAL(mpiNode.getHighCorner(i), 3);
	}

	bool coversAll[3] = {true};
	KDNode node(1, lowerCorner, higherCorner, 7, 0, coversAll, 2);
	node._deviationLowerBound = 2.0;
	node._deviation = 3.0;
	MPIKDNodePacked newMPINode = node.getMPIKDNode();
	ASSERT_EQUAL(newMPINode.getOwningProc(), mpiNode.getOwningProc());
	ASSERT_EQUAL(newMPINode.getNumProcs(), mpiNode.getNumProcs());
	ASSERT_EQUAL(newMPINode.getNodeID(), mpiNode.getNodeID());
	ASSERT_EQUAL(newMPINode.getFirstChildID(), mpiNode.getFirstChildID());
	ASSERT_EQUAL(newMPINode.getSecondChildID(), mpiNode.getSecondChildID());
	ASSERT_EQUAL(newMPINode.getLowCorner(0), mpiNode.getLowCorner(0));
	ASSERT_EQUAL(newMPINode.getLowCorner(1), mpiNode.getLowCorner(1));
	ASSERT_EQUAL(newMPINode.getLowCorner(2), mpiNode.getLowCorner(2));
	ASSERT_EQUAL(newMPINode.getHighCorner(0), mpiNode.getHighCorner(0));
	ASSERT_EQUAL(newMPINode.getHighCorner(1), mpiNode.getHighCorner(1));
	ASSERT_EQUAL(newMPINode.getHighCorner(2), mpiNode.getHighCorner(2));

	ASSERT_EQUAL(newMPINode.getDeviationLowerBound(), mpiNode.getDeviationLowerBound());
	ASSERT_EQUAL(newMPINode.getDeviation(), mpiNode.getDeviation());
	ASSERT_EQUAL(newMPINode.getLevel(), mpiNode.getLevel());
}

void KDNodeTest::testserializeDeserialize() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "KDNodeTest::testserializeDeserialize()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}
	int lowerEnd[] = {0, 0, 0};
	int upperEnd[] = {7, 7, 7};
	bool coversAll[] = {true, true, true};

	KDNode root(4, lowerEnd, upperEnd, 0, 0, coversAll, 0);
	root.buildKDTree();
//	std::cout << "root node:" << std::endl;
//	root.printTree();
	std::string name("KDNodeTest_serialize.out");
	root.serialize(name);

	KDNode result;
	result.deserialize(name);
//	std::cout << "Result: " << std::endl;
//	result.printTree();
	ASSERT_TRUE(root.equals(result));
}
