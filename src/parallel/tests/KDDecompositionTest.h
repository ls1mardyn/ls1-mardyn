/*
 * KDDecompositionTest.h
 *
 * @Date: 01.03.2012
 * @Author: eckhardw
 */

#ifndef KDDECOMPOSITIONTEST_H_
#define KDDECOMPOSITIONTEST_H_

#include "utils/TestWithSimulationSetup.h"
#include "parallel/KDNode.h"

#include <vector>

class KDDecompositionTest : public utils::TestWithSimulationSetup {

	TEST_SUITE(KDDecompositionTest);
	TEST_METHOD(testHaloCorrect);
	TEST_METHOD(testNoDuplicatedParticles);
	TEST_METHOD(testNoDuplicatedParticles2);
	TEST_METHOD(testNoLostParticles);
	TEST_METHOD(testNoLostParticles2);
	TEST_METHOD(testCompleteTreeInfo);
	TEST_METHOD(testRebalancingDeadlocks);
	TEST_METHOD(testbalanceAndExchange);
	TEST_SUITE_END();

public:

	KDDecompositionTest();

	virtual ~KDDecompositionTest();

	void testHaloCorrect();

	void testNoDuplicatedParticles();
	void testNoDuplicatedParticles2();

	void testNoLostParticles();
	void testNoLostParticles2();

	void testCompleteTreeInfo();

	/**
	 * Initial implementation of completeTreeInfo(). Kept to test the new / current
	 * implementation against it.
	 */
	void completeTreeInfo(KDNode*& root, KDNode*& ownArea, int ownRank);

	void testRebalancingDeadlocks();

	void testbalanceAndExchange();

private:

	void testNoDuplicatedParticlesFilename(const char * filename, double cutoff, double domainLength);

	void testNoLostParticlesFilename(const char * filename, double cutoff, double domainLength);

	/**
	 * init some random distribution
	 * @param v the _numParticlesPerCell entry of the KDD
	 * @param len the dimensions of the _numParticlesPerCell array of the KDD
	 */
	void setNumParticlesPerCell(std::vector<unsigned int> &v, int len[3]) const;

	/**
	 * set the entries to zero
	 * @param v the _numParticlesPerCell entry of the KDD
	 * @param total_len the total length of the _numParticlesPerCell array of the KDD
	 */
	void clearNumParticlesPerCell(std::vector<unsigned int> &v, int totalLen) const;

	unsigned f(double x, double y, double z, int N[3], const std::vector<double>& c) const;

	void Write_VTK_Structured_Points(unsigned *A, int N[3], std::string& filename) const;

	void initCoeffs(std::vector<double>& c) const;
	double myRand(double min, double max) const;

	std::vector<double> _currentCoeffs;
	std::vector<double> _oldCoeffs;
	int _rank;


};

#endif /* KDDECOMPOSITIONTEST_H_ */
