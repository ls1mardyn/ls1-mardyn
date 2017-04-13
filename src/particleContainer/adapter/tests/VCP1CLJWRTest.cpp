/*
 * VCP1CLJWRTest.cpp
 *
 *  Created on: 13 Apr 2017
 *      Author: tchipevn
 */

#include "VCP1CLJWRTest.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/VCP1CLJWR.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"

#ifdef MARDYN_WR
TEST_SUITE_REGISTRATION(VCP1CLJWRTest);
#else
#warning "The unit test for the WR force calculation is not executed in non-WR mode."
#endif


VCP1CLJWRTest::VCP1CLJWRTest() {
	test_log->info() << "Testing VCP1CLJWR cell processor against: " <<
#if VCP_VEC_TYPE==VCP_NOVEC
	"VectorizedCellProcessor with no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	"VectorizedCellProcessor with SSE intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	"VectorizedCellProcessor with AVX intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	"VectorizedCellProcessor with AVX2 intrinsics." << std::endl;
#endif
}

VCP1CLJWRTest::~VCP1CLJWRTest() {
}

void VCP1CLJWRTest::testForcePotentialCalculationU0() {
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}

	// U (r_ij) = 4 epsilon * ( (sigma / r_ij)^12 - (sigma / r_ij)^6 )
	// F (r_ij) = 24 epsilon * 1/r ( (sigma / r_ij)^6 - 2 * (sigma / r_ij)^12 )

	double forces[4][3] = { { -24, -24, 0 },
	                        {  24, -24, 0 },
	                        { -24,  24, 0 },
	                        {  24,  24, 0 }};

	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "ForceCalculationTestU0.inp", 1.1);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-8);

	VCP1CLJ_WR cellProcessor( *_domain, 1.1, 1.1);
	cellProcessor.setDtInv2m(1.0);
	container->traverseCells(cellProcessor);

	for (ParticleIterator m = container->iteratorBegin(); m != container->iteratorEnd(); ++m) {
		m->calcFM();
	}

	for (ParticleIterator m = container->iteratorBegin(); m != container->iteratorEnd(); ++m) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m->id() << " index i="<< i << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), forces[m->id()-1][i], m->F(i), 1e-8);
		}
	}

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(96, _domain->getLocalVirial(), 1e-8);

	delete container;
}

void VCP1CLJWRTest::testForcePotentialCalculationF0() {
}

void VCP1CLJWRTest::testLennardJonesVectorization() {
}
