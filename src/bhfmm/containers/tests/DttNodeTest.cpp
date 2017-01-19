/*
 * DttNodeTest.cpp
 *
 *  Created on: Nov 19, 2015
 *      Author: tchipevn
 */

#include "bhfmm/containers/tests/DttNodeTest.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"

TEST_SUITE_REGISTRATION(DttNodeTest);

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

	ParticleCell root;

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
  
	dtt::DttNode dummy(root, 0,ctr,globalDomainLength,globalDomainLength,orderOfExpansions,depth);
  
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
	ParticleContainer * container = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", 1.0);

	ParticleCell root;

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

	std::vector<double> shift;
	for(int i=0; i<3; i++){
		shift.push_back(4.78*(i+5));
	}
  for(int i=0; i<10000; i++){
//		root.convertAoSToSoACharge();
//		root.convertSoAToAoSCharge();
	}

	std::vector<Molecule*>& p_after = root.getParticlePointers();
	Molecule& m = *p_after[0]; 
	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong x coordinate",p_pos[0],m.r(0), 1e-12);
	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong y coordinate",p_pos[1],m.r(1), 1e-12);
	ASSERT_DOUBLES_EQUAL_MSG("P(1) wrong z coordinate",p_pos[2],m.r(2), 1e-12);
}

void DttNodeTest::testDepth(double cutoffRadius){
	double globalDomainLength[3] = {8., 8., 8.};
	double ctr[3] = {4., 4., 4.};
	int orderOfExpansions = 2;

	ParticleContainer * container = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", 1.0);

	ParticleCell root;

	Molecule * it;
	for(it = container->begin(); it != container->end(); it = container->next()) {
		root.addParticle(it);
	}
	
	int depth = log2((globalDomainLength[0] / cutoffRadius));
  
	dtt::DttNode dummy(root, 0,ctr,globalDomainLength,orderOfExpansions,depth);

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
	ParticleContainer * container = initializeFromFile(ParticleContainerFactory::LinkedCell, "FMMCharge.inp", 1.0);

	ParticleCell root;

	Molecule * it;
	for(it = container->begin(); it != container->end(); it = container->next()) {
		root.addParticle(it);
	}

	std::vector<ParticleCell> children(8, ParticleCell());

	int dummyOrder = 10;
	dtt::DttNode dummy(dummyOrder);
	dummy._ctr[0] = 4.;
	dummy._ctr[1] = 4.;
	dummy._ctr[2] = 4.;

	dummy.divideParticles(root,children); // TODO: root should be passed by reference?

	ASSERT_EQUAL_MSG("cell(0) should contain no particles", children[0].getMoleculeCount(), 0);
	ASSERT_EQUAL_MSG("cell(1) should contain one  particle", children[1].getMoleculeCount(), 1);
	ASSERT_EQUAL_MSG("cell(2) should contain no  particles", children[2].getMoleculeCount(), 0);
	ASSERT_EQUAL_MSG("cell(3) should contain no  particles", children[3].getMoleculeCount(), 0);
	ASSERT_EQUAL_MSG("cell(4) should contain no  particles", children[4].getMoleculeCount(), 0);
	ASSERT_EQUAL_MSG("cell(5) should contain no  particles", children[5].getMoleculeCount(), 0);
	ASSERT_EQUAL_MSG("cell(6) should contain two  particles", children[6].getMoleculeCount(), 2);
	ASSERT_EQUAL_MSG("cell(7) should contain one particles", children[7].getMoleculeCount(), 0);

	delete container;
}
