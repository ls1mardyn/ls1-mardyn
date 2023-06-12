/*
 * DttNodeTest.cpp
 *
 *  Created on: Nov 19, 2015
 *      Author: tchipevn
 */

#include "parallel/DomainDecompBase.h"
#include "bhfmm/containers/tests/DttNodeTest.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "bhfmm/containers/ParticleCellPointers.h"

#ifndef ENABLE_REDUCED_MEMORY_MODE
TEST_SUITE_REGISTRATION(DttNodeTest);
#else
#pragma message "Compilation info: DttNodeTest disabled in reduced memory mode"
#endif

DttNodeTest::DttNodeTest() {
	// TODO Auto-generated constructor stub
}

DttNodeTest::~DttNodeTest() {
	// TODO Auto-generated destructor stub
}

void DttNodeTest::testUpwardDownwardWithNoInteraction(){
/*	double globalDomainLength[3] = {8., 8., 8.};
	double ctr[3] = {4., 4., 4.};
	int orderOfExpansions = 2;

	ParticleContainer * container = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", 1.0);

	ParticleCellPointers root;

	double p_pos[3];
	bool check = true;

	Molecule * it;
	for(it = container->begin(); it != container->end(); it = container->next()) {
		root.addParticle(it);

		if(check){
			p_pos[0] = it->r(0);
			p_pos[1] = it->r(1);
			p_pos[2] = it->r(2);
		check = false;
		}
	}
	int depth = log2((globalDomainLength[0] / 1.0));

	bhfmm::DttNode dummy(root, 0,ctr,globalDomainLength,globalDomainLength,orderOfExpansions,depth);

	dummy.upwardPass();
//	dummy._mp_cell.local = dummy._mp_cell.multipole;
	dummy.downwardPass();

	dummy.getLeafParticleCell();

	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong x coordinate",p_pos[0],m.r(0), 1e-12);
	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong y coordinate",p_pos[1],m.r(1), 1e-12);
	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong z coordinate",p_pos[2],m.r(2), 1e-12);
*/
}

void DttNodeTest::testSoAConvertions(){
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testSoAConvertions()" << " not executed (rerun with only 1 Process!)"
				<< std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}
	ParticleContainer * container = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", 1.0);

	bhfmm::ParticleCellPointers root;
	double l[3], u[3];
	for (int d = 0; d < 3; ++d) {
		l[d] = container->getBoundingBoxMin(d);
		u[d] = container->getBoundingBoxMax(d);
	}
	root.setBoxMin(l);
	root.setBoxMax(u);

	double p_pos[3];
	bool check = true;

	for (auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
		root.addParticle(&(*it));

		if (check) {
			p_pos[0] = it->r(0);
			p_pos[1] = it->r(1);
			p_pos[2] = it->r(2);
			check = false;
		}
	}

	std::vector<double> shift;
	for (int i = 0; i < 3; i++) {
		shift.push_back(4.78 * (i + 5));
	}
	for (int i = 0; i < 10000; i++) {
//		root.convertAoSToSoACharge();
//		root.convertSoAToAoSCharge();
	}

	Molecule& m = root.moleculesAt(0);
	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong x coordinate",p_pos[0],m.r(0), 1e-12);
	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong y coordinate",p_pos[1],m.r(1), 1e-12);
	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong z coordinate",p_pos[2],m.r(2), 1e-12);
	delete container;
}

void DttNodeTest::testDepth(double cutoffRadius){
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testDepth()" << " not executed (rerun with only 1 Process!)"
				<< std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}
	double globalDomainLength[3] = {8., 8., 8.};
	double ctr[3] = {4., 4., 4.};
	bhfmm::Vector3<double> gDL_vec3(globalDomainLength);
	bhfmm::Vector3<double> ctr_vec3(ctr);

	int orderOfExpansions = 2;

	ParticleContainer * container = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", 1.0);

	std::vector<Molecule *> particles;

	for(auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
		particles.push_back(&(*it));
	}

	int depth = log2((globalDomainLength[0] / cutoffRadius));

	bhfmm::DttNode dummy(particles, 0,ctr,globalDomainLength,orderOfExpansions,depth);

	ASSERT_EQUAL_MSG("WRONG DEPTH!",dummy.getMaxDepth(),depth);

	delete container;
}

void DttNodeTest::testDepthAtRadius1(){
	testDepth(1.0);
}
void DttNodeTest::testDepthAtRadius2(){
	testDepth(2.0);
}
void DttNodeTest::testDepthAtRadius4(){
	testDepth(4.0);
}

void DttNodeTest::testDivideParticles() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "not executing testDivideParticles for more than 1 proc" << std::endl;
		return;
	}

	ParticleContainer * container = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", 1.0);

	std::vector<Molecule *> particles;

	for(auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
		particles.push_back(&(*it));
	}

	std::array<std::vector<Molecule*>, 8> children;

	int dummyOrder = 10;
	bhfmm::DttNode dummy(dummyOrder);
	dummy._ctr[0] = 4.;
	dummy._ctr[1] = 4.;
	dummy._ctr[2] = 4.;

	dummy.divideParticles(particles,children); // TODO: root should be passed by reference?

	ASSERT_EQUAL_MSG("cell(0) should contain one particle ", children[0].size(), 2ul);
	ASSERT_EQUAL_MSG("cell(1) should contain no  particle ", children[1].size(), 0ul);
	ASSERT_EQUAL_MSG("cell(2) should contain no  particles", children[2].size(), 0ul);
	ASSERT_EQUAL_MSG("cell(3) should contain no  particles", children[3].size(), 0ul);
	ASSERT_EQUAL_MSG("cell(4) should contain no  particles", children[4].size(), 0ul);
	ASSERT_EQUAL_MSG("cell(5) should contain no  particles", children[5].size(), 0ul);
	ASSERT_EQUAL_MSG("cell(6) should contain no  particles", children[6].size(), 0ul);
	ASSERT_EQUAL_MSG("cell(7) should contain two particles", children[7].size(), 1ul);

	delete container;
}
