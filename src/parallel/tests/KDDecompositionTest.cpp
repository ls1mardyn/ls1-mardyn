/*
 * KDDecompositionTest.cpp
 *
 * @Date: 01.03.2012
 * @Author: eckhardw
 */

#include "KDDecompositionTest.h"
#include "Domain.h"
#include "particleContainer/LinkedCells.h"

#include <sstream>
#include <cmath>

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
		_rank = _domainDecomposition->getRank();
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
		int fullSearchThreshold = 2;

		_domain->setGlobalLength(0, boxL);
		_domain->setGlobalLength(1, boxL);
		_domain->setGlobalLength(2, boxL);
		kdd = new KDDecomposition(cutOff, _domain, 1, fullSearchThreshold);

		double cellsInCutoffRadius = 1.0;
		double bBoxMin[3];
		double bBoxMax[3];
		for (int i = 0; i < 3; i++) {
			bBoxMin[i] = kdd->getBoundingBoxMin(i, _domain);
			bBoxMax[i] = kdd->getBoundingBoxMax(i, _domain);
		}
		moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, cutOff, cutOff, cellsInCutoffRadius);
		moleculeContainer->update();
		kdd->_steps = 0;
		_rank = kdd->_rank;


		srand(42);
	}

	// TEST

	const int numReps = 10;
	if (_rank == 0) {
		cout << "running " << numReps << " repetitions" << std::endl;
	}
	for (int i = 0; i < numReps; ++i) {

		// initialize currentCoeffs
		initCoeffs(_currentCoeffs);

		KDNode * newDecompRoot = NULL;
		KDNode * newOwnLeaf = NULL;

		setNumParticlesPerCell(kdd->_numParticlesPerCell, kdd->_globalCellsPerDim);

		kdd->constructNewTree(newDecompRoot, newOwnLeaf, moleculeContainer);

		clearNumParticlesPerCell(kdd->_numParticlesPerCell, kdd->_globalNumCells);

		kdd->barrier();
		bool isOK =
				kdd->migrateParticles(*newDecompRoot, *newOwnLeaf, moleculeContainer);


		if(not isOK and _rank == 0)
		{
			std::stringstream fOld, fNew;

			fOld << "old_" << i << ".vtu";
			kdd->_decompTree->plotNode(fOld.str());

			fNew << "new_" << i << ".vtu";
			newDecompRoot->plotNode(fNew.str());

			cout << "current coeffs: " << std::endl;
			cout << setprecision(17);
			for (unsigned int j = 0; j < _currentCoeffs.size(); ++j) {
				cout << _currentCoeffs[j] << std::endl;
			}
			cout << std::endl;
			cout << "old coeffs: " << std::endl;
			for (unsigned int j = 0; j < _oldCoeffs.size(); ++j) {
				cout << _oldCoeffs[j] << std::endl;
			}
		}

		kdd->barrier();

		ASSERT_TRUE_MSG("Deadlock!", isOK);
		if (not isOK)
			MPI_Abort(MPI_COMM_WORLD, 1);

		delete kdd->_decompTree;
		kdd->_decompTree = newDecompRoot;
		kdd->_ownArea = newOwnLeaf;

		_oldCoeffs = _currentCoeffs;
		_currentCoeffs.clear();

		if(_rank == 0) {
			cout << "repetition " << i << " completed" << std::endl;
		}
	}


	// SHUTDOWN
	delete moleculeContainer;
	delete kdd;

}

void KDDecompositionTest::initCoeffs(std::vector<double>& c) const {
	for(int i = 0; i < 10; ++i)
		c.push_back(myRand(-1.0, 1.0));
	for(int i = 10; i < 22; ++i) {
		c.push_back(myRand(-2.0, 2.0));
	}
}

double KDDecompositionTest::myRand(double min, double max) const {
	double ret;
	double x =((double)rand()/(double)RAND_MAX); // in [0,1]
	ret = min + (max - min) * x;

	return ret;
}

void KDDecompositionTest::setNumParticlesPerCell(unsigned int * v, int len[3]) const {

	for (int z = 0; z < len[2]; ++z) {
		for (int y = 0; y < len[1]; ++y) {
			for (int x = 0; x < len[0]; ++x) {
				v[(z * len[1] + y) * len[0] + x] = f(x, y, z, len, _currentCoeffs);
			}
		}
	}

//	static int i = 0;
//	char fname[100];
//	sprintf(fname, "partsInCell_%i.vtk", i);
//	if(_rank == 0) {
//		Write_VTK_Structured_Points(v, len, fname);
//	}
//	++i;
}

unsigned KDDecompositionTest::f(double x, double y, double z, int N[3], const std::vector<double>& c) const {
	unsigned ret;

	double s =
//			c[0] * x + c[1] * y + c[2] * z +
//			c[3] * x * y + c[4] * y * z + c[5] * x * z +
//			c[6] * x * y * z +
//			c[7] * x / y + c[8] * y / z + c[9] * z / x +
			c[10] * sin(c[16] * x * M_PI / N[0] + c[0]) +
			c[11] * cos(c[17] * x * M_PI / N[0] + c[1]) +
			c[12] * sin(c[18] * y * M_PI / N[1] + c[2]) +
			c[13] * cos(c[19] * y * M_PI / N[1] + c[3]) +
			c[14] * sin(c[20] * z * M_PI / N[2] + c[4]) +
			c[15] * cos(c[21] * z * M_PI / N[2] + c[5]);
	s *= (c[6] + 20);
	s += 20;

	ret = fabs(s);
	ret %= 1000;

	return ret;
}

void KDDecompositionTest::clearNumParticlesPerCell(unsigned int * v, int totalLen) const {
	for (int i = 0; i < totalLen; ++i)
		v[i] = 0;
}

void KDDecompositionTest::Write_VTK_Structured_Points(unsigned *A, int N[3], const char *filename) const {
	FILE *out;
	out = fopen(filename, "w");
	fprintf(out, "# vtk DataFile Version 3.0\n");
	fprintf(out, "Laplace\n");
	fprintf(out, "ASCII\n");
	fprintf(out, "DATASET STRUCTURED_POINTS\n");
	fprintf(out, "DIMENSIONS %d %d %d\n", N[0], N[1], N[2]);
	fprintf(out, "ORIGIN %d %d %d\n", 0, 0, 0);
	fprintf(out, "SPACING %d %d %d\n", 1, 1, 1);
	fprintf(out, "POINT_DATA %d\n", N[0] * N[1] * N[2]);
	fprintf(out, "SCALARS nMols int 1\n");
	fprintf(out, "LOOKUP_TABLE default\n");
	for (int i = 0; i < N[0] * N[1] * N[2]; i++) {
		fprintf(out, "%i\n", A[i]);
	}
	fclose(out);
}
