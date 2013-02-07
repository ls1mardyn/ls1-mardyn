///*
// * KDDecompositionTest.cpp
// *
// * @Date: 01.03.2012
// * @Author: eckhardw
// */
//
//#include "KDDecompositionTest.h"
//#include "Domain.h"
//
//TEST_SUITE_REGISTRATION(KDDecompositionTest);
//
//using namespace std;
//
//KDDecompositionTest::KDDecompositionTest() {
//}
//
//KDDecompositionTest::~KDDecompositionTest() {
//}
//
//void KDDecompositionTest::testCompleteTreeInfo() {
//
//	if (_domainDecomposition->getNumProcs() < 9) {
//		_domain->setGlobalLength(0, 50);
//		_domain->setGlobalLength(1, 50);
//		_domain->setGlobalLength(2, 50);
//
//		/*
//		 * create two times the same tree and communicate the first, which should
//		 * not change it.
//		 */
//		int lowerEnd[] = {0, 0, 0};
//		int upperEnd[] = {3, 3, 3};
//		bool coversAll[] = {true, true, true};
//		KDNode* root = new KDNode(_domainDecomposition->getNumProcs(), lowerEnd, upperEnd, 0, 0, coversAll, 0);
//		root->buildKDTree();
//		KDNode* ownArea = root->findAreaForProcess(_rank);
//
//		KDNode result(_domainDecomposition->getNumProcs(), lowerEnd, upperEnd, 0, 0, coversAll, 0);
//		result.buildKDTree();
//
//		KDDecomposition decomposition(1.0, _domain, 1.0, 10);
//		decomposition.completeTreeInfo(root, ownArea);
//		ASSERT_TRUE(result.equals(*root));
//
//	} else {
//		Log::global_log->warning() << "KDDecompositionTest::testCompleteTreeInfo():"
//				"only executed with 8 or less procs!"<< std::endl;
//	}
//}

