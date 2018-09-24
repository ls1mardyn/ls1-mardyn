/**
 * @file AutoPasContainerTest.h
 * @author seckler
 * @date 19.09.18
 */

#pragma once

#include "ParticleContainerTest.h"

class AutoPasContainerTest : public ParticleContainerTest {

	TEST_SUITE(AutoPasContainerTest);
	TEST_METHOD(testConstructor);
/*
    TEST_METHOD(testInsertion);
	TEST_METHOD(testMoleculeIteration);
	TEST_METHOD(testUpdateAndDeleteOuterParticles);
	TEST_METHOD(testUpdateAndDeleteOuterParticlesH2O);
	TEST_METHOD(testUpdateAndDeleteOuterParticles8Particles);
	TEST_METHOD(testMoleculeBeginNextEndDeleteCurrent);
	TEST_METHOD(testTraversalMethods);

	TEST_METHOD(testCellBorderAndFlagManager);
 */

	TEST_SUITE_END();

public:
	AutoPasContainerTest(){};
	void testConstructor();
};




