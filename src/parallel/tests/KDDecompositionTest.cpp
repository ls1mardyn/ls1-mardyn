/*
 * KDDecompositionTest.cpp
 *
 * @Date: 01.03.2012
 * @Author: eckhardw
 */

#include "KDDecompositionTest.h"
#include "Domain.h"
#include "particleContainer/LinkedCells.h"

TEST_SUITE_REGISTRATION(KDDecompositionTest);

using namespace std;

KDDecompositionTest::KDDecompositionTest() {
}

KDDecompositionTest::~KDDecompositionTest() {
}

void KDDecompositionTest::testCompleteTreeInfo() {

	if (_domainDecomposition->getNumProcs() < 9) {
		_domain->setGlobalLength(0, 50);
		_domain->setGlobalLength(1, 50);
		_domain->setGlobalLength(2, 50);

		/*
		 * create two times the same tree and communicate the first, which should
		 * not change it.
		 */
		int lowerEnd[] = {0, 0, 0};
		int upperEnd[] = {3, 3, 3};
		bool coversAll[] = {true, true, true};
		KDNode* root = new KDNode(_domainDecomposition->getNumProcs(), lowerEnd, upperEnd, 0, 0, coversAll, 0);
		root->buildKDTree();
		KDNode* ownArea = root->findAreaForProcess(_rank);

		KDNode result(_domainDecomposition->getNumProcs(), lowerEnd, upperEnd, 0, 0, coversAll, 0);
		result.buildKDTree();

		KDDecomposition decomposition(1.0, _domain, 1.0, 10);
		decomposition.completeTreeInfo(root, ownArea);
		ASSERT_TRUE(result.equals(*root));

	} else {
		Log::global_log->warning() << "KDDecompositionTest::testCompleteTreeInfo():"
				"only executed with 8 or less procs!"<< std::endl;
	}
}

void KDDecompositionTest::testRebalancingDeadlocks() {

	// INIT
	KDDecomposition * kdd;
	LinkedCells * moleculeContainer;
	{
		const double boxL = 1241.26574;
		const double cutOff = 26.4562;

		_domain->setGlobalLength(0, boxL);
		_domain->setGlobalLength(1, boxL);
		_domain->setGlobalLength(2, boxL);
		kdd = new KDDecomposition(cutOff, _domain, 1, 3);

		int cellsInCutoffRadius;
		double bBoxMin[3];
		double bBoxMax[3];
		for (int i = 0; i < 3; i++) {
			bBoxMin[i] = kdd->getBoundingBoxMin(i, _domain);
			bBoxMax[i] = kdd->getBoundingBoxMax(i, _domain);
		}
		moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, cutOff, cutOff, 1.0);
		moleculeContainer->update();
		kdd->_steps = 0;
	}

	// TEST

	KDNode * newDecompRoot = NULL;
	KDNode * newOwnLeaf = NULL;

	setNumParticlesPerCell(kdd->_numParticlesPerCell, kdd->_globalCellsPerDim);
	kdd->constructNewTree(newDecompRoot, newOwnLeaf);

	clearNumParticlesPerCell(kdd->_numParticlesPerCell, kdd->_globalNumCells);
	kdd->migrateParticles(*newDecompRoot, *newOwnLeaf, moleculeContainer);

	delete kdd->_decompTree;
	kdd->_decompTree = newDecompRoot;
	kdd->_ownArea = newOwnLeaf;



	// SHUTDOWN
	delete moleculeContainer;
	delete kdd;

}

void KDDecompositionTest::setNumParticlesPerCell(unsigned int * v, int len[3]) const {
	// dummy implementation test
	for (int z = 0; z < len[2]; ++z) {
		for (int y = 0; y < len[1]; ++y) {
			for (int x = 0; x < len[0]; ++x) {
				v[(z * len[1] + y) * len[0] + x] = x + y + z;
			}
		}
	}
}

void KDDecompositionTest::clearNumParticlesPerCell(unsigned int * v, int totalLen) const {
	for (int i = 0; i < totalLen; ++i)
		v[i] = 0;
}
