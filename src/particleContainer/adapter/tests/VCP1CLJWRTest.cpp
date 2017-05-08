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
#pragma message "Compilation info: The unit test for the WR force calculation is not executed in non-WR mode."
#endif


VCP1CLJWRTest::VCP1CLJWRTest() {
#ifdef MARDYN_WR
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
#endif /* MARDYN_WR */
}

VCP1CLJWRTest::~VCP1CLJWRTest() {
}

void VCP1CLJWRTest::testForcePotentialCalculationU0() {
#ifdef MARDYN_WR
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "VCP1CLJWRTest::testForcePotentialCalculationU0()"
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
	cellProcessor.setDtInvm(1.0);
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
#endif /* MARDYN_WR */
}

void VCP1CLJWRTest::testForcePotentialCalculationF0() {
#ifdef MARDYN_WR
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "VCP1CLJWRTest::testForcePotentialCalculationF0()"
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
	cellProcessor.setDtInvm(1.0);
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
#endif /* MARDYN_WR */
}

// free function
void VCP1CLJWRTest__initFullCellSoA(const ParticleCell_WR & cell_wr, CellDataSoA& fullSoA) {
#ifdef MARDYN_WR

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

		//for better readability:
		constexpr ConcatenatedSites<vcp_real_calc>::SiteType LJC = ConcatenatedSites<vcp_real_calc>::SiteType::LJC;
		typedef ConcatenatedSites<vcp_real_calc>::CoordinateType Coordinate;
		typedef CellDataSoA::QuantityType QuantityType;

		fullSoA.getBegin(QuantityType::MOL_POSITION, LJC, Coordinate::X)[ind] = M.r(0);
		fullSoA.getBegin(QuantityType::MOL_POSITION, LJC, Coordinate::Y)[ind] = M.r(1);
		fullSoA.getBegin(QuantityType::MOL_POSITION, LJC, Coordinate::Z)[ind] = M.r(2);
		fullSoA.getBegin(QuantityType::CENTER_POSITION, LJC, Coordinate::X)[ind] = M.r(0);
		fullSoA.getBegin(QuantityType::CENTER_POSITION, LJC, Coordinate::Y)[ind] = M.r(1);
		fullSoA.getBegin(QuantityType::CENTER_POSITION, LJC, Coordinate::Z)[ind] = M.r(2);
		fullSoA._ljc_id[ind] = 0;

		// clear FM
		std::array<vcp_real_calc, 3> clearance = { 0., 0., 0. };
		fullSoA.setTriplet(clearance, QuantityType::FORCE, LJC, ind);
		fullSoA.setTriplet(clearance, QuantityType::VIRIAL, LJC, ind);	}
#endif /* MARDYN_WR */
}

void VCP1CLJWRTest::testProcessCell() {
#ifdef MARDYN_WR
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "VCP1CLJWRTest::testProcessCell()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}

	double ScenarioCutoff = 35.0;
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "VectorizationLennardJones1CLJ.inp", ScenarioCutoff);
	for (ParticleIterator m = container->iteratorBegin(); m != container->iteratorEnd(); ++m) {
		for (int d = 0; d < 3; ++d) {
			m->setv(d, 0.0);
		}
	}

	LinkedCells * linkedCells = dynamic_cast<LinkedCells*>(container);

	VCP1CLJ_WR vcp_WR( *_domain,  ScenarioCutoff, ScenarioCutoff);
	vcp_WR.setDtInvm(1.0);
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

	size_t numMolecules = cell_wr.getMoleculeCount();
	for (size_t i = 0; i < numMolecules; ++i) {
		double WR_f_x = cell_wr.moleculesAt(i).F(0);
		double WR_f_y = cell_wr.moleculesAt(i).F(1);
		double WR_f_z = cell_wr.moleculesAt(i).F(2);

		std::array<vcp_real_calc, 3> triple = full_SoA.getTriplet(CellDataSoA::QuantityType::FORCE, ConcatenatedSites<vcp_real_calc>::SiteType::LJC, i); //LJC equals the beginning of data
		double full_f_x = static_cast<double>triple[0];
		double full_f_y = static_cast<double>triple[1];
		double full_f_z = static_cast<double>triple[2];
		ASSERT_DOUBLES_EQUAL_MSG("force x should have been equal.", full_f_x, WR_f_x, 1.0e-20);
		ASSERT_DOUBLES_EQUAL_MSG("force y should have been equal.", full_f_y, WR_f_y, 1.0e-20);
		ASSERT_DOUBLES_EQUAL_MSG("force z should have been equal.", full_f_z, WR_f_z, 1.0e-20);
	}

	delete container;
#endif /* MARDYN_WR */
}

void VCP1CLJWRTest::testProcessCellPair() {
#ifdef MARDYN_WR
	if (_domainDecomposition->getNumProcs() != 1) {
		test_log->info() << "VCP1CLJWRTest::testProcessCellPair()"
				<< " not executed (rerun with only 1 Process!)" << std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
		return;
	}

	// copy-paste cause I'm lazy and have no particular time for unit tests.
	double ScenarioCutoff = 35.0;
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "VectorizationLennardJones1CLJ.inp", ScenarioCutoff);
	for (ParticleIterator m = container->iteratorBegin(); m != container->iteratorEnd(); ++m) {
		for (int d = 0; d < 3; ++d) {
			m->setv(d, 0.0);
		}
	}

	LinkedCells * linkedCells = dynamic_cast<LinkedCells*>(container);

	VCP1CLJ_WR vcp_WR( *_domain,  ScenarioCutoff, ScenarioCutoff);
	vcp_WR.setDtInvm(1.0);
	VectorizedCellProcessor vcp_full(*_domain,  ScenarioCutoff, ScenarioCutoff);

	// get an inner cell
	double innerPoint[3] = {0.1, 0.1, 0.1};
	unsigned long firstCellIndex = linkedCells->getCellIndexOfPoint(innerPoint);
	ParticleCell_WR& cell_wr1 = linkedCells->getCell(firstCellIndex);
	ParticleCell_WR& cell_wr2 = linkedCells->getCell(firstCellIndex + 1);


	CellDataSoA full_SoA1(0,0,0,0,0);
	CellDataSoA full_SoA2(0,0,0,0,0);
	VCP1CLJWRTest__initFullCellSoA(cell_wr1, full_SoA1);
	VCP1CLJWRTest__initFullCellSoA(cell_wr2, full_SoA2);

	vcp_WR.initTraversal();
	vcp_WR.processCellPair(cell_wr1, cell_wr2);
	vcp_WR.endTraversal();


	double WR_Upot = _domain->getLocalUpot();
	double WR_Virial = _domain->getLocalVirial();

	vcp_full.initTraversal();
	const bool CalculateMacroscopic = true;
	const bool ApplyCutoff = true;
	vcp_full._calculatePairs<CellPairPolicy_<ApplyCutoff>, CalculateMacroscopic, MaskGatherC>(full_SoA2, full_SoA1);
	vcp_full.endTraversal();

	double full_Upot = _domain->getLocalUpot();
	double full_Virial = _domain->getLocalVirial();

	ASSERT_DOUBLES_EQUAL(full_Upot, WR_Upot, 1.0e-17);
	ASSERT_DOUBLES_EQUAL(full_Virial, WR_Virial, 1.0e-17);

	size_t numMolecules1 = cell_wr1.getMoleculeCount();
	for (size_t i = 0; i < numMolecules1; ++i) {
		double WR_f_x = cell_wr1.moleculesAt(i).F(0);
		double WR_f_y = cell_wr1.moleculesAt(i).F(1);
		double WR_f_z = cell_wr1.moleculesAt(i).F(2);

		std::array<vcp_real_calc, 3> triple = full_SoA1.getTriplet(CellDataSoA::QuantityType::FORCE, ConcatenatedSites<vcp_real_calc>::SiteType::LJC, i); //LJC equals the beginning of data
		double full_f_x = static_cast<double>(triple[0]);
		double full_f_y = static_cast<double>(triple[1])
		double full_f_z = static_cast<double>(triple[2])
		ASSERT_DOUBLES_EQUAL_MSG("force x should have been equal.", full_f_x, WR_f_x, 1.0e-20);
		ASSERT_DOUBLES_EQUAL_MSG("force y should have been equal.", full_f_y, WR_f_y, 1.0e-20);
		ASSERT_DOUBLES_EQUAL_MSG("force z should have been equal.", full_f_z, WR_f_z, 1.0e-20);
	}

	size_t numMolecules2 = cell_wr2.getMoleculeCount();
	for (size_t i = 0; i < numMolecules2; ++i) {
		double WR_f_x = cell_wr2.moleculesAt(i).F(0);
		double WR_f_y = cell_wr2.moleculesAt(i).F(1);
		double WR_f_z = cell_wr2.moleculesAt(i).F(2);

		std::array<vcp_real_calc, 3> triple = full_SoA2.getTriplet(CellDataSoA::QuantityType::FORCE, ConcatenatedSites<vcp_real_calc>::SiteType::LJC, i); //LJC equals the beginning of data
		double full_f_x = static_cast<double>triple[0];
		double full_f_y = static_cast<double>triple[1];
		double full_f_z = static_cast<double>triple[2];
		ASSERT_DOUBLES_EQUAL_MSG("force x should have been equal.", full_f_x, WR_f_x, 1.0e-20);
		ASSERT_DOUBLES_EQUAL_MSG("force y should have been equal.", full_f_y, WR_f_y, 1.0e-20);
		ASSERT_DOUBLES_EQUAL_MSG("force z should have been equal.", full_f_z, WR_f_z, 1.0e-20);
	}

	delete container;
#endif /* MARDYN_WR */
}
