/*
 * ParticleCellRMMTest.cpp
 *
 *  Created on: 12 Oct 2017
 *      Author: tchipevn
 */

#include "ParticleCellRMMTest.h"
#include "particleContainer/ParticleCellRMM.h"

TEST_SUITE_REGISTRATION(ParticleCellRMMTest);

ParticleCellRMMTest::ParticleCellRMMTest() {
	// TODO Auto-generated constructor stub

}

ParticleCellRMMTest::~ParticleCellRMMTest() {
	// TODO Auto-generated destructor stub
}

void ParticleCellRMMTest::testSizeOfIs64() {
	std::string msg = "sizeof(ParticleCellRMM) should be 64 bytes. Please do not add member fields to this class or any of the base classes for the moment. If really necessary, try to make this work or contact the TUM developers.";
	ASSERT_EQUAL_MSG(msg, sizeof(ParticleCellRMM), 64ul);
}
