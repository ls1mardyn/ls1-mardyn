/*
 * VCP1CLJWRTest.cpp
 *
 *  Created on: 13 Apr 2017
 *      Author: tchipevn
 */

#include "VCP1CLJWRTest.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/adapter/VCP1CLJWR.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/adapter/vectorization/MaskGatherChooser.h"

#ifdef MARDYN_WR
TEST_SUITE_REGISTRATION(VCP1CLJWRTest);
#else
#warning "The unit test for the WR force calculation is not executed in non-WR mode."
#endif


VCP1CLJWRTest::VCP1CLJWRTest() {
	test_log->info() << "Testing VCP1CLJWR cell processor against: "
#if VCP_VEC_TYPE==VCP_NOVEC
	<< "VectorizedCellProcessor with no intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
	<< "VectorizedCellProcessor with SSE intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	<< "VectorizedCellProcessor with AVX intrinsics." << std::endl;
#elif VCP_VEC_TYPE==VCP_VEC_AVX2
	<< "VectorizedCellProcessor with AVX2 intrinsics." << std::endl;
#endif
	;
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
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "DomainDecompositionTest::testExchangeMolecules1Proc()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}


	// U (r_ij) = 4 epsilon * ( (sigma / r_ij)^12 - (sigma / r_ij)^6 )
	// F (r_ij) = 24 epsilon * 1/r ( (sigma / r_ij)^6 - 2 * (sigma / r_ij)^12 )
	// -> potential U = -1 per particle pair.
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "ForceCalculationTestF0.inp", 1.3);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-8);

	VCP1CLJ_WR cellProcessor( *_domain,  1.3, 1.3);
	cellProcessor.setDtInv2m(1.0);
	container->traverseCells(cellProcessor);

	for (ParticleIterator m = container->iteratorBegin(); m != container->iteratorEnd(); ++m) {
		m->calcFM();
	}

	for (ParticleIterator m = container->iteratorBegin(); m != container->iteratorEnd(); ++m) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m->id() << " index i="<< i << " F[i]=" << m->F(i) << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), 0.0, m->F(i), 1e-7);
		}
	}

	ASSERT_DOUBLES_EQUAL(-4, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-6);

	delete container;
}

// free function
void VCP1CLJWRTest__initFullCellSoA(const ParticleCell_WR & cell_wr, CellDataSoA& fullSoA) {

//	for (int i = 0; i < numMols; ++i) {
////		FullMolecule m = FullMolecule(cell_wr.moleculesAtConst(i));
//		TODO: create an std vector of FullMolecules and SoA for two cells
//		TODO: compute forces via _calculatePairs
//	}

	// Determine the total number of centers.
	size_t numMolecules = cell_wr.getMoleculeCount();
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	size_t nDipoles = 0;
	size_t nQuadrupoles = 0;

	for (size_t m = 0;  m < numMolecules; ++m) {
		const Molecule_WR& Mol = cell_wr.moleculesAtConst(m);
		nLJCenters += Mol.numLJcenters();
		nCharges += Mol.numCharges();
		nDipoles += Mol.numDipoles();
		nQuadrupoles += Mol.numQuadrupoles();
	}
	mardyn_assert(nCharges == 0);
	mardyn_assert(nDipoles == 0);
	mardyn_assert(nQuadrupoles == 0);

	// Construct the SoA.
	fullSoA.resize(numMolecules,nLJCenters,nCharges,nDipoles,nQuadrupoles);

	size_t iLJCenters = 0;
	size_t iCharges = 0;
	size_t iDipoles = 0;
	size_t iQuadrupoles = 0;

	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < numMolecules; ++i) {
		const Molecule & M = cell_wr.moleculesAtConst(i);
		const size_t mol_ljc_num = M.numLJcenters();
		const size_t mol_charges_num = M.numCharges();
		const size_t mol_dipoles_num = M.numDipoles();
		const size_t mol_quadrupoles_num = M.numQuadrupoles();

		fullSoA._mol_ljc_num[i] = mol_ljc_num;
		fullSoA._mol_charges_num[i] = mol_charges_num;
		fullSoA._mol_dipoles_num[i] = mol_dipoles_num;
		fullSoA._mol_quadrupoles_num[i] = mol_quadrupoles_num;

		fullSoA._mol_pos.x(i) = M.r(0);
		fullSoA._mol_pos.y(i) = M.r(1);
		fullSoA._mol_pos.z(i) = M.r(2);

		mardyn_assert( M.numLJcenters() == 1);
		const unsigned ind = i;
		fullSoA.ljc_m_r_xBegin()[ind] = M.r(0);
		fullSoA.ljc_m_r_yBegin()[ind] = M.r(1);
		fullSoA.ljc_m_r_zBegin()[ind] = M.r(2);
		fullSoA.ljc_r_xBegin()[ind] = M.r(0);
		fullSoA.ljc_r_yBegin()[ind] = M.r(1);
		fullSoA.ljc_r_zBegin()[ind] = M.r(2);
		fullSoA._ljc_id[ind] = 0;

		// clear FM
		fullSoA.ljc_f_xBegin()[ind] = M.r(0);
		fullSoA.ljc_f_yBegin()[ind] = M.r(1);
		fullSoA.ljc_f_zBegin()[ind] = M.r(2);
		fullSoA.ljc_V_xBegin()[ind] = M.r(0);
		fullSoA.ljc_V_yBegin()[ind] = M.r(1);
		fullSoA.ljc_V_zBegin()[ind] = M.r(2);
	}
}

void VCP1CLJWRTest::testLennardJonesVectorization() {
	double ScenarioCutoff = 35.0;
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "VectorizationLennardJones1CLJ.inp", ScenarioCutoff);
	LinkedCells * linkedCells = dynamic_cast<LinkedCells*>(container);

	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-8);

	VCP1CLJ_WR vcp_WR( *_domain,  ScenarioCutoff, ScenarioCutoff);
	vcp_WR.setDtInv2m(1.0);
	VectorizedCellProcessor vcp_full(*_domain,  ScenarioCutoff, ScenarioCutoff);

	// get an inner cell
	double innerPoint[3] = {0.1, 0.1, 0.1};
	unsigned long firstCellIndex = linkedCells->getCellIndexOfPoint(innerPoint);
	ParticleCell_WR& cell_wr = linkedCells->getCell(firstCellIndex);
	CellDataSoA full_SoA(0,0,0,0,0);
	VCP1CLJWRTest__initFullCellSoA(cell_wr, full_SoA);

	vcp_WR.initTraversal();
	vcp_WR.processCell(cell_wr);
	vcp_WR.endTraversal();

	double WR_Upot = _domain->getLocalUpot();
	double WR_Virial = _domain->getLocalVirial();

	vcp_full.initTraversal();
	const bool CalculateMacroscopic = true;
	const bool ApplyCutoff = true;
	vcp_full._calculatePairs<SingleCellPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(full_SoA, full_SoA);
	vcp_full.endTraversal();

	double full_Upot = _domain->getLocalUpot();
	double full_Virial = _domain->getLocalVirial();

	ASSERT_DOUBLES_EQUAL(full_Upot, WR_Upot, 1.0e-17);
	ASSERT_DOUBLES_EQUAL(full_Virial, WR_Virial, 1.0e-17);
	test_log->info() << "WR computed Upot :" << WR_Upot << std::endl;
	test_log->info() << "WR computed Virial :" << WR_Virial << std::endl;
	test_log->info() << "num Molecules " << cell_wr.getMoleculeCount() << std::endl;


#if 0
	for (ParticleIterator m = container->iteratorBegin(); m != container->iteratorEnd(); ++m) {
		for (int i = 0; i < 3; i++) {
			std::stringstream str;
			str << "Molecule id=" << m->id() << " index i="<< i << " F[i]=" << m->F(i) << std::endl;
			ASSERT_DOUBLES_EQUAL_MSG(str.str(), 0.0, m->F(i), 1e-7);
		}
	}

	ASSERT_DOUBLES_EQUAL(-4, _domain->getLocalUpot(), 1e-8);
	ASSERT_DOUBLES_EQUAL(0.0, _domain->getLocalVirial(), 1e-6);
#endif

	delete container;

}
