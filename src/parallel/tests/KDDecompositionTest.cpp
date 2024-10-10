/*
 * KDDecompositionTest.cpp
 *
 * @Date: 01.03.2012
 * @Author: eckhardw
 */

#include "KDDecompositionTest.h"
#include "Domain.h"
#ifdef MARDYN_AUTOPAS
#include "particleContainer/AutoPasContainer.h"
#endif
#include "particleContainer/LinkedCells.h"
#include "io/ASCIIReader.h"
#include "parallel/NeighbourCommunicationScheme.h"
#include <utils/mardyn_assert.h>

#include <sstream>
#include <cmath>
#include <string>
#include <fstream>

TEST_SUITE_REGISTRATION(KDDecompositionTest);


KDDecompositionTest::KDDecompositionTest() :
		_rank(0) {
}

KDDecompositionTest::~KDDecompositionTest() {
}

void KDDecompositionTest::testNoDuplicatedParticlesFilename(
		const char * filename, double cutoff, double domainLength) {
	KDDecomposition * kdd;
	// original pointer will be deleted by tearDown()
	_domain->setGlobalLength(0, domainLength);
	_domain->setGlobalLength(1, domainLength);
	_domain->setGlobalLength(2, domainLength);
	kdd = new KDDecomposition(cutoff, 3, 1, 4);
	kdd->init(_domain);
	_domainDecomposition = kdd;
	_rank = kdd->_rank;
	ParticleContainer* container = initializeFromFile(
			ParticleContainerFactory::LinkedCell, filename, cutoff);

	kdd->initCommunicationPartners(cutoff, _domain, container);

	int numMols = container->getNumberOfParticles();
	kdd->collCommInit(1);
	kdd->collCommAppendInt(numMols);
	kdd->collCommAllreduceSum();
	numMols = kdd->collCommGetInt();
	kdd->collCommFinalize();

	_domainDecomposition->balanceAndExchange(0., true, container, _domain);
	// will rebalance, we thus need a reduce
	container->deleteOuterParticles();

	int newNumMols = container->getNumberOfParticles();

	kdd->collCommInit(1);
	kdd->collCommAppendInt(newNumMols);
	kdd->collCommAllreduceSum();
	newNumMols = kdd->collCommGetInt();
	kdd->collCommFinalize();

	ASSERT_EQUAL(numMols, newNumMols);

	delete _domainDecomposition;
	delete container;
}

/**
 * checks if the halo is correct (not duplicate and no lost)
 */
void KDDecompositionTest::testHaloCorrect() {
	KDDecomposition * kdd;
	// original pointer will be deleted by tearDown()
	double globalLength = 13.;
	_domain->setGlobalLength(0, globalLength);
	_domain->setGlobalLength(1, globalLength);
	_domain->setGlobalLength(2, globalLength);
	double cutoff = 2.5;
	kdd = new KDDecomposition(cutoff, 1, 5, 4);
	kdd->init(_domain);
	_domainDecomposition = kdd;
	_rank = kdd->_rank;

	// this file contains a regular particle distribution with particles at the positions (i,k,l)
	// with i,k,l in {1., 2., ... , 11., 12.}
	ParticleContainer* container = initializeFromFile(
			ParticleContainerFactory::LinkedCell, "haloCorrect_regular.inp",
			cutoff);

	for (int d = 0; d < 3; ++d) {
		kdd->_neighbourCommunicationScheme->setCoverWholeDomain(d, kdd->_ownArea->_coversWholeDomain[d]);
	}
	kdd->initCommunicationPartners(cutoff, _domain, container);
	_domainDecomposition->balanceAndExchange(0., true, container, _domain);
	// this should give us halos at -1., -2., 14., 15.
	// we will thus iterate over

	double pos[3];
	double boxMin[3] = { container->getBoundingBoxMin(0),
			container->getBoundingBoxMin(1), container->getBoundingBoxMin(2) };
	double boxMax[3] = { container->getBoundingBoxMax(0),
			container->getBoundingBoxMax(1), container->getBoundingBoxMax(2) };
	for (pos[0] = -2.; pos[0] < 15.5; pos[0] += 1.) {
		if (pos[0] < boxMin[0] - cutoff or pos[0] > boxMax[0] + cutoff
				or pos[0] == 0 or pos[0] == 13) {
			//not in our box or pos == 0. or 13.
			continue;
		}
		for (pos[1] = -2.; pos[1] < 15.5; pos[1] += 1.) {
			if (pos[1] < boxMin[1] - cutoff or pos[1] > boxMax[1] + cutoff
					or pos[1] == 0 or pos[1] == 13) {
				//not in our box or pos == 0. or 13.
				continue;
			}
			for (pos[2] = -2.; pos[2] < 15.5; pos[2] += 1.) {
				if (pos[2] < boxMin[2] - cutoff or pos[2] > boxMax[2] + cutoff
						or pos[2] == 0 or pos[2] == 13) {
					//not in our box or pos == 0 or 13.
					continue;
				}
				if ((pos[0] >= boxMin[0] and pos[0] < boxMax[0])
						and (pos[1] >= boxMin[1] and pos[1] < boxMax[1])
						and (pos[2] >= boxMin[2] and pos[2] < boxMax[2])) {
					// inside, so normal molecule
					continue;
				}
				auto iter = container->getMoleculeAtPosition(pos);
				std::visit([](auto iter) { ASSERT_TRUE_MSG("halo molecule not present", iter.isValid()); }, iter);
			}
		}
	}

	delete _domainDecomposition;
	delete container;
}

void KDDecompositionTest::testNoDuplicatedParticles() {
	testNoDuplicatedParticlesFilename("H20_NaBr_0.01_T_293.15_DD.inp", 5.0,
			58.5389);
}
void KDDecompositionTest::testNoDuplicatedParticles2() {
	testNoDuplicatedParticlesFilename("H20_NaBr_0.01_T_293.15_DD_2.inp", 5.0,
			58.5);
}

void KDDecompositionTest::testNoLostParticlesFilename(const char * filename,
		double cutoff, double domainLength) {
	// original pointer will be deleted by tearDown()
	KDDecomposition * kdd;
	// original pointer will be deleted by tearDown()
	_domain->setGlobalLength(0, domainLength);
	_domain->setGlobalLength(1, domainLength);
	_domain->setGlobalLength(2, domainLength);
	kdd = new KDDecomposition(cutoff, 3, 1, 4);
	kdd->init(_domain);
	_domainDecomposition = kdd;
	_rank = kdd->_rank;

	ParticleContainer* container = initializeFromFile(
			ParticleContainerFactory::LinkedCell, filename, cutoff);

	for (int d = 0; d < 3; ++d) {
		kdd->_neighbourCommunicationScheme->setCoverWholeDomain(d, kdd->_ownArea->_coversWholeDomain[d]);
	}
	kdd->initCommunicationPartners(cutoff, _domain, container);

	int numMols = container->getNumberOfParticles();
	kdd->collCommInit(1);
	kdd->collCommAppendInt(numMols);
	kdd->collCommAllreduceSum();
	numMols = kdd->collCommGetInt();
	kdd->collCommFinalize();

	double bBoxMin[3];
	double bBoxMax[3];
	for (int dim = 0; dim < 3; dim++) {
		bBoxMin[dim] = 0.;
		bBoxMax[dim] = _domain->getGlobalLength(dim);
	}
	std::set<unsigned long> lower[3]; // the id of particles that were close to the lower boundary in the specific dimension are stored here
	std::set<unsigned long> upper[3]; // the id of particles that were close to the upper boundary in the specific dimension are stored here

#if defined(_OPENMP)
#pragma omp parallel
#endif
	{

		std::set<unsigned long> lower_thread[3]; // the id of particles that were close to the lower boundary in the specific dimension are stored here
		std::set<unsigned long> upper_thread[3]; // the id of particles that were close to the upper boundary in the specific dimension are stored here

		for (auto m = container->iterator(ParticleIterator::ALL_CELLS); m.isValid();
				++m) {
			for (int dim = 0; dim < 3; dim++) {
				if (m->r(dim) < bBoxMin[dim] + cutoff * 0.5) {
					// we shift particles close to the lower boundary to outside of the lower boundary.
					// in this case they are put to the smallest (in abs values) negative representable number
					// i.e. 2^(-149) = -1.4013e-45 for float resp. 4.94066e-324 for double
					m->setr(dim,
							std::nexttoward((vcp_real_calc) bBoxMin[dim],
									bBoxMin[dim] - 1.f));
					lower_thread[dim].insert(m->getID());
				}
				if (m->r(dim) > bBoxMax[dim] - cutoff * 0.5) {
					// We shift particles close to the upper boundary to outside of the upper boundary.
					// In this case they are put at minimum to boundingBoxMax, as this is no longer inside of the domain.
					// If the float representation of the maximum is less than the double representation, the next bigger floating point representation is used.
					// Otherwise the maximum is used.
					vcp_real_calc r =
							(float) bBoxMax[dim] >= bBoxMax[dim] ?
									bBoxMax[dim] :
									std::nexttoward(
											(vcp_real_calc) bBoxMax[dim],
											bBoxMax[dim] + 1.f);
					m->setr(dim, r);
					upper_thread[dim].insert(m->getID());
				}
			}
		}

#if defined(_OPENMP)
#pragma omp critical
#endif
		{
			for (int d = 0; d < 3; ++d) {
				for (auto it = lower_thread[d].begin();
						it != lower_thread[d].end(); ++it)
					lower[d].insert(*it);

				for (auto it = upper_thread[d].begin();
						it != upper_thread[d].end(); ++it)
					upper[d].insert(*it);
			}
		}
	}

	container->update();

	for(auto iter = container->iterator(ParticleIterator::Type::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter){
		bool found = false;
		auto id = iter->getID();
		for(unsigned short dim = 0; dim < 3; ++dim){
			if(lower[dim].count(id) or upper[dim].count(id)){
				found = true;
			}
		}
		ASSERT_EQUAL_MSG("shifted molecule still in inner cells", found, false);
	}

	//needed to properly exchange the particles. In the first step leaving particles are normally not exchanged...
	dynamic_cast<KDDecomposition*>(_domainDecomposition)->_steps++;

	_domainDecomposition->balanceAndExchange(0., true, container, _domain);
	container->deleteOuterParticles();

	int newNumMols = container->getNumberOfParticles();

	kdd->collCommInit(1);
	kdd->collCommAppendInt(newNumMols);
	kdd->collCommAllreduceSum();
	newNumMols = kdd->collCommGetInt();
	kdd->collCommFinalize();

	//_domain->writeCheckpoint("dump.txt", container, _domainDecomposition, false);
	ASSERT_EQUAL(numMols, newNumMols);

	for (auto m = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); m.isValid(); ++m) {
		for (int dim = 0; dim < 3; dim++) {
			if (lower[dim].count(m->getID())) {
				// We make sure, that these particles are now at the top part of the domain.
				ASSERT_TRUE(m->r(dim) >= bBoxMax[dim] - cutoff / 2.);
			} else if (upper[dim].count(m->getID())) {
				// We make sure, that these particles are now at the lower part of the domain.
				ASSERT_TRUE(m->r(dim) <= bBoxMin[dim] + cutoff / 2.);
			}
		}
	}

	delete _domainDecomposition;
	delete container;
}

void KDDecompositionTest::testNoLostParticles() {
	testNoLostParticlesFilename("H20_NaBr_0.01_T_293.15_DD.inp", 3.0, 58.5389);
}
void KDDecompositionTest::testNoLostParticles2() {
	testNoLostParticlesFilename("H20_NaBr_0.01_T_293.15_DD_2.inp", 3.0, 58.5);
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
		int lowerEnd[] = { 0, 0, 0 };
		int upperEnd[] = { 3, 3, 3 };
		bool coversAll[] = { true, true, true };
		KDNode* root = new KDNode(_domainDecomposition->getNumProcs(), lowerEnd,
				upperEnd, 0, 0, coversAll, 0);
		root->buildKDTree();
		_rank = _domainDecomposition->getRank();
		KDNode* ownArea = root->findAreaForProcess(_rank);

		KDNode result(_domainDecomposition->getNumProcs(), lowerEnd, upperEnd,
				0, 0, coversAll, 0);
		result.buildKDTree();

		KDDecomposition decomposition(1.0, 1, 1.0, 10);
		decomposition.init(_domain);
		KDNode * toCleanUp = root;
		decomposition.completeTreeInfo(root, ownArea);
		delete toCleanUp;
		ASSERT_TRUE(result.equals(*root));
		delete root;

	} else {
		Log::global_log->warning()
				<< "KDDecompositionTest::testCompleteTreeInfo():"
						"only executed with 8 or less procs!" << std::endl;
	}
}

void KDDecompositionTest::testRebalancingDeadlocks() {

	// INIT
	KDDecomposition * kdd;
	ParticleContainer * moleculeContainer;
	{
		const double boxL = 1241.26574;
		const double cutOff = 26.4562;
		int fullSearchThreshold = 2;

		_domain->setGlobalLength(0, boxL);
		_domain->setGlobalLength(1, boxL);
		_domain->setGlobalLength(2, boxL);
		kdd = new KDDecomposition(cutOff, 1, 1, fullSearchThreshold);
		kdd->init(_domain);

		double bBoxMin[3];
		double bBoxMax[3];
		for (int i = 0; i < 3; i++) {
			bBoxMin[i] = kdd->getBoundingBoxMin(i, _domain);
			bBoxMax[i] = kdd->getBoundingBoxMax(i, _domain);
		}
#ifndef MARDYN_AUTOPAS
		moleculeContainer = new LinkedCells(bBoxMin, bBoxMax, cutOff);
#else
		moleculeContainer = new AutoPasContainer(cutOff);
		moleculeContainer->rebuild(bBoxMin, bBoxMax);
#endif
		moleculeContainer->update();
		kdd->_steps = 0;
		_rank = kdd->_rank;

		srand(42);
	}

	// TEST

	const int numReps = 10;
	if (_rank == 0) {
		//cout << "running " << numReps << " repetitions" << std::endl;
	}
	for (int i = 0; i < numReps; ++i) {

		// initialize currentCoeffs
		initCoeffs(_currentCoeffs);

		KDNode * newDecompRoot = nullptr;
		KDNode * newOwnLeaf = nullptr;

		setNumParticlesPerCell(kdd->_numParticlesPerCell,
				kdd->_globalCellsPerDim);

		kdd->constructNewTree(newDecompRoot, newOwnLeaf, moleculeContainer);

		clearNumParticlesPerCell(kdd->_numParticlesPerCell,
				kdd->_globalNumCells);

		kdd->barrier();
		bool isOK = kdd->migrateParticles(*newDecompRoot, *newOwnLeaf,
				moleculeContainer, _domain);

		if (not isOK and _rank == 0) {
			std::stringstream fOld, fNew;

			fOld << "old_" << i << ".vtu";
			kdd->_decompTree->plotNode(fOld.str());

			fNew << "new_" << i << ".vtu";
			newDecompRoot->plotNode(fNew.str());

			std::cout << "current coeffs: " << std::endl;
			std::cout << std::setprecision(17);
			for (double _currentCoeff : _currentCoeffs) {
				std::cout << _currentCoeff << std::endl;
			}
			std::cout << std::endl;
			std::cout << "old coeffs: " << std::endl;
			for (double _oldCoeff : _oldCoeffs) {
				std::cout << _oldCoeff << std::endl;
			}
		}

		kdd->barrier();

		ASSERT_TRUE_MSG("Deadlock!", isOK);
		if (not isOK) {
			std::ostringstream error_message;
			error_message << "[KDDecompositionTest] Deadlock detected." << std::endl;
			MARDYN_EXIT(error_message);
		}
		delete kdd->_decompTree;
		kdd->_decompTree = newDecompRoot;
		kdd->_ownArea = newOwnLeaf;

		_oldCoeffs = _currentCoeffs;
		_currentCoeffs.clear();

		if (_rank == 0) {
			//cout << "repetition " << i << " completed" << std::endl;
		}
	}

	// SHUTDOWN
	delete moleculeContainer;
	delete kdd;

}

void KDDecompositionTest::testbalanceAndExchange() {

	// INIT
	KDDecomposition * kdd;

	const double cutOff = 3.5;
	int fullSearchThreshold = 2;

	ASCIIReader inputReader;
	std::string fileName2 = getTestDataFilename("DomainDecompBase.inp");
	inputReader.setPhaseSpaceHeaderFile(fileName2.c_str());
	inputReader.setPhaseSpaceFile(fileName2.c_str());
	inputReader.readPhaseSpaceHeader(_domain, 1.0);

	kdd = new KDDecomposition(cutOff, 1, 1, fullSearchThreshold);
	kdd->init(_domain);
	_domainDecomposition = kdd;

	_rank = kdd->_rank;

	ParticleContainer* moleculeContainer = initializeFromFile(
			ParticleContainerFactory::LinkedCell, "DomainDecompBase.inp",
			cutOff);

	kdd->initCommunicationPartners(cutOff, _domain, moleculeContainer);

	moleculeContainer->update();

	// TEST

	const int numReps = 10;
	if (_rank == 0) {
		//cout << "running " << numReps << " repetitions" << std::endl;
	}
	for (int i = 0; i < numReps; ++i) {
		kdd->balanceAndExchange(1.0, true, moleculeContainer, _domain);
		moleculeContainer->updateMoleculeCaches();
	}
	delete moleculeContainer;
	delete kdd;

	// SHUTDOWN

}

void KDDecompositionTest::initCoeffs(std::vector<double>& c) const {
	for (int i = 0; i < 10; ++i)
		c.push_back(myRand(-1.0, 1.0));
	for (int i = 10; i < 22; ++i) {
		c.push_back(myRand(-2.0, 2.0));
	}
}

double KDDecompositionTest::myRand(double min, double max) const {
	double ret;
	double x = ((double) rand() / (double) RAND_MAX); // in [0,1]
	ret = min + (max - min) * x;

	return ret;
}

void KDDecompositionTest::setNumParticlesPerCell(std::vector<unsigned int> &v,
		int len[3]) const {

	for (int z = 0; z < len[2]; ++z) {
		for (int y = 0; y < len[1]; ++y) {
			for (int x = 0; x < len[0]; ++x) {
				v[(z * len[1] + y) * len[0] + x] = f(x, y, z, len,
						_currentCoeffs);
			}
		}
	}

//	if(_rank == 0) {
//		std::string fname = "partsInCell_0.vtk";
//		Write_VTK_Structured_Points(v, len, fname);
//	}
}

unsigned KDDecompositionTest::f(double x, double y, double z, int N[3],
		const std::vector<double>& c) const {
	unsigned ret;

	double s =
//			c[0] * x + c[1] * y + c[2] * z +
//			c[3] * x * y + c[4] * y * z + c[5] * x * z +
//			c[6] * x * y * z +
//			c[7] * x / y + c[8] * y / z + c[9] * z / x +
			c[10] * sin(c[16] * x * M_PI / N[0] + c[0])
					+ c[11] * cos(c[17] * x * M_PI / N[0] + c[1])
					+ c[12] * sin(c[18] * y * M_PI / N[1] + c[2])
					+ c[13] * cos(c[19] * y * M_PI / N[1] + c[3])
					+ c[14] * sin(c[20] * z * M_PI / N[2] + c[4])
					+ c[15] * cos(c[21] * z * M_PI / N[2] + c[5]);
	s *= (c[6] + 20);
	s += 20;

	ret = fabs(s);
	ret %= 1000;

	return ret;
}

void KDDecompositionTest::clearNumParticlesPerCell(std::vector<unsigned int> &v,
		int totalLen) const {
	for (int i = 0; i < totalLen; ++i)
		v[i] = 0;
}

void KDDecompositionTest::Write_VTK_Structured_Points(unsigned *A, int N[3], std::string& filename) const {
	std::ofstream ofs(filename, std::ios::out);
	ofs << "# vtk DataFile Version 3.0" << std::endl;
	ofs << "Laplace" << std::endl;
	ofs << "ASCII" << std::endl;
	ofs << "DATASET STRUCTURED_POINTS" << std::endl;
	ofs << "DIMENSIONS " << N[0] << " " << N[1] << " " << N[2] << std::endl;
	ofs << "ORIGIN 0 0 0"<< std::endl;
	ofs << "SPACING 1 1 1" << std::endl;
	ofs << "POINT_DATA " << N[0]*N[1]*N[2] << std::endl;
	ofs << "SCALARS nMols int 1" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < N[0] * N[1] * N[2]; i++) {
			ofs << A[i] << std::endl;
	}
	ofs.close();
}
