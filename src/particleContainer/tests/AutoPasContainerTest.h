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
	TEST_METHOD(testUpdate);

	TEST_SUITE_END();

public:
	AutoPasContainerTest() = default;
	void testConstructor();
	void testUpdate();
};




