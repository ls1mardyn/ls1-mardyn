/*
 * LinkedCellsTest.cpp
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

class LinkedCellsTest;
#include "LinkedCellsTest.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif
#include "particleContainer/adapter/CellProcessor.h"
#include <vector>

#include "ensemble/PressureGradient.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/LinkedCellTraversals/HalfShellTraversal.h"
#include "particleContainer/TraversalTuner.h"

TEST_SUITE_REGISTRATION(LinkedCellsTest);

LinkedCellsTest::LinkedCellsTest() {

}

LinkedCellsTest::~LinkedCellsTest() {
}

void LinkedCellsTest::testUpdateAndDeleteOuterParticlesFilename(const char * filename, double cutoff) {
	// original pointer will be deleted by tearDown()
	_domainDecomposition = new DomainDecompBase();

	LinkedCells* container = static_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
			filename, cutoff));
	int numMols = container->getNumberOfParticles();

	_domainDecomposition->exchangeMolecules(container, _domain);
	container->deleteOuterParticles();

	int newNumMols = container->getNumberOfParticles();
//	_domain->writeCheckpoint("dump.txt", container, _domainDecomposition);
	ASSERT_EQUAL(numMols, newNumMols);

	delete _domainDecomposition;
	delete container;
}

void LinkedCellsTest::testUpdateAndDeleteOuterParticlesH2O() {
	testUpdateAndDeleteOuterParticlesFilename("H20_NaBr_0.01_T_293.15.inp", 27.0);
}

void LinkedCellsTest::testUpdateAndDeleteOuterParticles8Particles() {
	testUpdateAndDeleteOuterParticlesFilename("LinkedCells.inp", 1.0);
}

void LinkedCellsTest::testMoleculeBeginNextEndDeleteCurrent() {
#ifdef ENABLE_REDUCED_MEMORY_MODE
	global_log->warning() << "LinkedCellsTest::testMoleculeBeginNextEndDeleteCurrent() needs to be redone in REDUCED_MEMORY_MODE (it is currently disabled)."
	<< std::endl;
	//TODO: take a real scenario and proper LinkedCells object, don't be so lazy.
	return;
#endif
	// NOTE: we do not open an OpenMP parallel region!
	// Hence, this test is always executed sequentially!

	Molecule dummyMolecule1(1, &_components[0], 1.0, 1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule2(2, &_components[0], 2.0, 2.0, 2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule3(3, &_components[0], 3.0, 3.0, 3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule4(4, &_components[0], 5.1, 5.1, 5.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	LinkedCells LC;
	ParticleIterator molIt;

	std::vector < ParticleCell > &cells = LC._cells;

	cells.resize(5);

	// cell[0] is empty

	// cell[1] contains 1,2,3
	{
		double l[3] = { 0., 0., 0. }, u[3] = { 5., 5., 5. };
		cells[1].setBoxMin(l);
		cells[1].setBoxMax(u);
	}
	cells[1].addParticle(dummyMolecule1);
	cells[1].addParticle(dummyMolecule2);
	cells[1].addParticle(dummyMolecule3);

	// cell[2] is empty

	// cell[3] contains 4
	{
		double l[3] = { 5., 5., 5. }, u[3] = { 5.5, 5.5, 5.5 };
		cells[3].setBoxMin(l);
		cells[3].setBoxMax(u);
	}
	cells[3].addParticle(dummyMolecule4);

	// cell[4] is empty

	// BEGIN:
	molIt = LC.iteratorBegin();
	ASSERT_TRUE_MSG("begin()", molIt->id() == 1ul);
	ASSERT_TRUE_MSG("end()", molIt != LC.iteratorEnd());

	// NEXT:
	++molIt;
	ASSERT_TRUE_MSG("next() within cell", molIt->id() == 2ul);
	ASSERT_TRUE_MSG("end()", molIt != LC.iteratorEnd());
	++molIt;
	ASSERT_TRUE_MSG("next() within cell", molIt->id() == 3ul);
	ASSERT_TRUE_MSG("end()", molIt != LC.iteratorEnd());
	++molIt;
	ASSERT_TRUE_MSG("next() across cells", molIt->id() == 4ul);
	ASSERT_TRUE_MSG("end()", molIt != LC.iteratorEnd());
	++molIt;
	ASSERT_TRUE_MSG("next() arrive at end()", molIt == LC.iteratorEnd());

	// DELETECURRENT:
	molIt = LC.iteratorBegin();

	molIt.deleteCurrentParticle();
	++molIt;
	ASSERT_EQUAL_MSG("delete() within cell", 3ul, molIt->id()); // 3 copied in place of 1
	molIt.deleteCurrentParticle();
	++molIt;
	ASSERT_TRUE_MSG("delete() within cell", molIt->id() == 2ul); // 2 copied in place of 3
	molIt.deleteCurrentParticle();
	++molIt;
	ASSERT_TRUE_MSG("delete() across cells", molIt->id() == 4ul); // cell 1 became empty, we advanced to cell 3
	molIt.deleteCurrentParticle();
	++molIt;
	ASSERT_TRUE_MSG("delete() last", molIt == LC.iteratorEnd()); // cell 4 became empty, we arrived at end()
}

void LinkedCellsTest::testParticleIteratorBeginNextEndParticleIteratorSequential() {
#ifdef ENABLE_REDUCED_MEMORY_MODE
	global_log->warning() << "LinkedCellsTest::testParticleIteratorBeginNextEndParticleIteratorSequential() needs to be redone in REDUCED_MEMORY_MODE (it is currently disabled)."
	<< std::endl;
	//TODO: take a real scenario and proper LinkedCells object, don't be so lazy.
	return;
#endif
	// NOTE: we do not open an OpenMP parallel region!
	// Hence, this test is always executed sequentially!

	Molecule dummyMolecule1(1, &_components[0], 1.0, 1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule2(2, &_components[0], 2.0, 2.0, 2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule3(3, &_components[0], 3.0, 3.0, 3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule4(4, &_components[0], 5.1, 5.1, 5.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	LinkedCells LC;

	std::vector < ParticleCell > &cells = LC._cells;

	cells.resize(5);

	// cell[0] is empty

	// cell[1] contains 1,2,3
	{
		double l[3] = { 0., 0., 0. }, u[3] = { 5., 5., 5. };
		cells[1].setBoxMin(l);
		cells[1].setBoxMax(u);
	}
	cells[1].addParticle(dummyMolecule1);
	cells[1].addParticle(dummyMolecule2);
	cells[1].addParticle(dummyMolecule3);

	// cell[2] is empty

	// cell[3] contains 4
	{
		double l[3] = { 5., 5., 5. }, u[3] = { 5.5, 5.5, 5.5 };
		cells[3].setBoxMin(l);
		cells[3].setBoxMax(u);
	}
	cells[3].addParticle(dummyMolecule4);

	// cell[4] is empty

	ParticleIterator begin = LC.iteratorBegin();
	ParticleIterator end = LC.iteratorEnd();

	// test that molecule IDs are: 1, 2, 3, 4, in this order
	// and begin and end work correctly
	unsigned long uID = 1;
	for (ParticleIterator mol = begin; mol != end; ++mol, ++uID) {
		ASSERT_EQUAL(uID, mol->id());
	}
}

#if 0
void LinkedCellsTest::testGetHaloBoundaryParticlesDirection() {
#if 0
	moleculeContainer->getRegion(_regionLow, _regionHigh, particles);
	moleculeContainer->getHaloParticlesDirection(_signedDirection, part2, false);
	moleculeContainer->getBoundaryParticlesDirection(_signedDirection, part2);
#endif

	double rCut = 40.;
	// make sure we have a DomainDecompBase

	_domainDecomposition = new DomainDecompBase();
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, "Ethan_equilibrated.inp", rCut);

	double regionLow[3] = {-rCut, -rCut, -rCut};
	double regionHigh[3] = {rCut, 571.607759, 571.607759};

	std::vector<Molecule *> partsRegion;
	container->getRegion(regionLow, regionHigh, partsRegion);

	int signedDirection = -1;
	std::vector<Molecule *> partsHaloAndBoundary;
	container->getHaloParticlesDirection(signedDirection, partsHaloAndBoundary, false);
	container->getBoundaryParticlesDirection(signedDirection, partsHaloAndBoundary);

	std::vector<Molecule *> partsFiltered;

	for (int i = 0; i < partsHaloAndBoundary.size(); ++i) {
		Molecule * m = partsHaloAndBoundary[i];
		if (m->r(0) < rCut)
		partsFiltered.push_back(m);
	}

	ASSERT_EQUAL(partsRegion.size(), partsFiltered.size());
	delete _domainDecomposition;

}

#endif /* 0 */

class CellProcessorStub: public CellProcessor {
public:
	CellProcessorStub(size_t num_cells) :
			CellProcessor(0., 0.), _cellProcessCount(num_cells, 0), _cellPairProcessCount(num_cells), _num_cells(
					num_cells) {
	}
	virtual void initTraversal() {
	}

	virtual void preprocessCell(ParticleCell& /*cell*/) {
	}

	virtual void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false) { // does this need a bool
		_cellPairProcessCount[cell1.getCellIndex()][cell2.getCellIndex()] += sign;
		_cellPairProcessCount[cell2.getCellIndex()][cell1.getCellIndex()] += sign;  // newton 3
	}

	virtual void processCell(ParticleCell& cell) {
		_cellProcessCount[cell.getCellIndex()] += sign;
	}

	virtual double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) {
		return 0.;
	}

	virtual void postprocessCell(ParticleCell& /*cell*/) {
	}

	virtual void endTraversal() {

	}

	void checkZero() {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
		for (int i = 0; i < _num_cells; ++i) {
			ASSERT_EQUAL(_cellProcessCount[i], 0);

			std::map<unsigned long, int> & m = _cellPairProcessCount[i];
			for (auto& j : m) {
				ASSERT_EQUAL(j.second, 0);
			}
		} /* end of parallel region */
	}

	void checkOnlyInner(LinkedCells* container) {
#if defined(_OPENMP)
#pragma omp parallel for
#endif
		for (int i = 0; i < _num_cells; i++) {
			if (_cellProcessCount[i] != 0) {
				ASSERT_TRUE(container->getCellReference(i).isInnerCell());
			}

			for (auto& pair : _cellPairProcessCount[i]) {
				if (pair.second != 0) {
					ASSERT_TRUE(
							container->getCellReference(i).isInnerCell()
									&& container->getCellReference(pair.first).isInnerCell());
				}
			}
		}/* end of parallel region */
	}

	void inverseSign() {
		sign *= -1;
	}
private:
	std::vector<int> _cellProcessCount;
	std::vector<std::map<unsigned long, int>> _cellPairProcessCount;
	int _num_cells;
	int sign = 1;
};

void LinkedCellsTest::testTraversalMethods() {
	const char* filename = "VectorizationMultiComponentMultiPotentials.inp";
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 5.);
	int* boxWidthInNumCells = dynamic_cast<LinkedCells*>(container)->boxWidthInNumCells();
	int haloWidthInNumCells = container->getHaloWidthNumCells();
	size_t numCells = (boxWidthInNumCells[0] + 2 * haloWidthInNumCells)
			* (boxWidthInNumCells[1] + 2 * haloWidthInNumCells) * (boxWidthInNumCells[2] + 2 * haloWidthInNumCells);
	CellProcessorStub cpStub(numCells);

	container->traverseCells(cpStub);
	cpStub.inverseSign();
	container->traversePartialInnermostCells(cpStub, 0, 1);
	container->traverseNonInnermostCells(cpStub);
	cpStub.checkZero();
	cpStub.inverseSign();

	container->traverseCells(cpStub);
	cpStub.inverseSign();
	container->traversePartialInnermostCells(cpStub, 0, 3);
	container->traversePartialInnermostCells(cpStub, 1, 3);
	container->traversePartialInnermostCells(cpStub, 2, 3);
	container->traverseNonInnermostCells(cpStub);
	cpStub.checkZero();
	cpStub.inverseSign();

	container->traversePartialInnermostCells(cpStub, 0, 1);
	cpStub.checkOnlyInner(static_cast<LinkedCells*>(container));
	cpStub.inverseSign();
	container->traversePartialInnermostCells(cpStub, 0, 1);
	cpStub.inverseSign();

	container->traversePartialInnermostCells(cpStub, 0, 3);
	container->traversePartialInnermostCells(cpStub, 1, 3);
	container->traversePartialInnermostCells(cpStub, 2, 3);
	cpStub.checkOnlyInner(static_cast<LinkedCells*>(container));
	cpStub.inverseSign();
	container->traversePartialInnermostCells(cpStub, 0, 3);
	container->traversePartialInnermostCells(cpStub, 1, 3);
	container->traversePartialInnermostCells(cpStub, 2, 3);
	cpStub.inverseSign();
	delete container;
}

//void LinkedCellsTest::testHalfShell() {
//	//TODO: ___Extract to separate test class
//	//------------------------------------------------------------
//	// Setup
//	//------------------------------------------------------------
//
//	if (_domainDecomposition->getNumProcs() != 1) {
//		test_log->info() << "LinkedCellsTest::testHalfShell()"
//				<< " not executed (rerun with only 1 Process!)" << std::endl;
//		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
//		return;
//	}
//
//	auto domainDecomposition = new DomainDecompBase();
//	auto filename = "LinkedCellsHS.inp";
//	auto cutoff = 1;
//
//	LinkedCells* containerHS = dynamic_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
//			filename, cutoff));
//	containerHS->_traversalTuner->selectedTraversal = TraversalTuner<ParticleCell>::traversalNames::HS;
//	containerHS->_traversalTuner->findOptimalTraversal();
//	containerHS->initializeTraversal();
//
//	LinkedCells* container = dynamic_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
//			filename, cutoff));
//
//	auto vectorizedCellProcessor = new VectorizedCellProcessor(*_domain, cutoff, cutoff);
//
//
//	//------------------------------------------------------------
//	//  Calculate forces for FS and HS and compare
//	//------------------------------------------------------------
//
//	doHSTest(domainDecomposition, vectorizedCellProcessor, container, containerHS);
//
//	//------------------------------------------------------------
//	// Cleanup
//	//------------------------------------------------------------
//
//	delete domainDecomposition;
//	delete vectorizedCellProcessor;
//}
//
//void LinkedCellsTest::testMidpoint() {
//
//	//TODO: ___Extract to separate test class
//	//------------------------------------------------------------
//	// Setup
//	//------------------------------------------------------------
//
//	if (_domainDecomposition->getNumProcs() != 1) {
//		test_log->info() << "LinkedCellsTest::testMidpoint()"
//				<< " not executed (rerun with only 1 Process!)" << std::endl;
//		std::cout << "numProcs:" << _domainDecomposition->getNumProcs() << std::endl;
//		return;
//	}
//
//	auto domainDecomposition = new DomainDecompBase();
//	auto filename = "LinkedCellsHS.inp";
//	auto cutoff = 1;
//
//	LinkedCells* containerMP = dynamic_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
//			filename, cutoff));
//	containerMP->_traversalTuner->selectedTraversal = TraversalTuner<ParticleCell>::traversalNames::MP;
//	containerMP->_traversalTuner->findOptimalTraversal();
//	containerMP->initializeTraversal();
//
//	LinkedCells* container = dynamic_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
//			filename, cutoff));
//
//	auto vectorizedCellProcessor = new VectorizedCellProcessor(*_domain, cutoff, cutoff);
//
//
//	//------------------------------------------------------------
//	//  Calculate forces for FS and MP and compare
//	//------------------------------------------------------------
//
//	doHSTest(domainDecomposition, vectorizedCellProcessor, container, containerMP);
//
//	//------------------------------------------------------------
//	// Cleanup
//	//------------------------------------------------------------
//
//	delete domainDecomposition;
//	delete vectorizedCellProcessor;
//}

// new tests here

void LinkedCellsTest::testHalfShellMPIIndirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "indirect", "hs");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "indirect", "hs");
}

void LinkedCellsTest::testHalfShellMPIDirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "direct", "hs");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "direct", "hs");
}

void LinkedCellsTest::testMidpointMPIIndirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "indirect", "mp");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "indirect", "mp");
}

void LinkedCellsTest::testMidpointMPIDirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "direct", "mp");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "direct", "mp");
}

void LinkedCellsTest::doForceComparisonTest(std::string inputFile,
		TraversalTuner<ParticleCell>::traversalNames traversal, unsigned cellsInCutoff, std::string neighbourCommScheme,
		std::string commScheme) {

	//------------------------------------------------------------
	// Setup
	//------------------------------------------------------------

#ifdef ENABLE_MPI
	auto domainDecompositionFS = new DomainDecomposition();
	auto domainDecompositionTest = new DomainDecomposition();
	domainDecompositionFS->setCommunicationScheme(neighbourCommScheme, "fs");
	domainDecompositionTest->setCommunicationScheme(neighbourCommScheme, commScheme);

#else
	auto domainDecompositionFS = new DomainDecompBase();
	auto domainDecompositionTest = new DomainDecompBase();
#endif
	auto filename = inputFile.c_str();
	auto cutoff = 3.5;

	LinkedCells* containerTest = dynamic_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
			filename, cutoff));
	containerTest->_traversalTuner->selectedTraversal = traversal;
	containerTest->initializeTraversal();
	containerTest->_cellsInCutoff = cellsInCutoff;

	LinkedCells* container = dynamic_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
			filename, cutoff));

	// Legacy Processor for debugging
	// assert in potforce.h l. 483 has to be disabled

	/*
	 auto pphandler = new ParticlePairs2PotForceAdapter(*_domain);
	 auto pphandler2 = new ParticlePairs2PotForceAdapter(*_domain);
	 auto cellProc = new LegacyCellProcessor(cutoff, cutoff, pphandler);
	 auto cellProc2 = new LegacyCellProcessor(cutoff, cutoff, pphandler2);
	 */

	auto cellProc = new VectorizedCellProcessor(*_domain, cutoff, cutoff);
	auto cellProc2 = new VectorizedCellProcessor(*_domain, cutoff, cutoff);

#ifdef ENABLE_MPI
	domainDecompositionFS->initCommunicationPartners(cutoff, _domain);
	domainDecompositionTest->initCommunicationPartners(cutoff, _domain);
#endif

	//------------------------------------------------------------
	//  Calculate forces for FS and TestTraversal and compare
	//------------------------------------------------------------
	//------------------------------------------------------------
	// Prepare molecule containers
	//------------------------------------------------------------
	container->update();
	containerTest->update();
	bool forceRebalancing = false;
	domainDecompositionFS->balanceAndExchange(0.,forceRebalancing, container, _domain);
	domainDecompositionTest->balanceAndExchange(0.,forceRebalancing, containerTest, _domain);
	container->updateMoleculeCaches();
	containerTest->updateMoleculeCaches();
	//------------------------------------------------------------
	// Do calculation with FS
	//------------------------------------------------------------
	{
		container->traverseCells(*cellProc);
		// calculate forces
		const ParticleIterator begin = container->iteratorBegin();
		const ParticleIterator end = container->iteratorEnd();
		for (auto i = begin; i != end; ++i) {
			i->calcFM();
		}
	}
	//------------------------------------------------------------
	// Do calculation with TestTraversal
	//------------------------------------------------------------
	{
		containerTest->traverseCells(*cellProc2);
		// calculate forces
		const ParticleIterator& begin = containerTest->iteratorBegin();
		const ParticleIterator& end = containerTest->iteratorEnd();
		for (ParticleIterator i = begin; i != end; ++i) {
			i->calcFM();
//			std::cout << "r: " << i->r(0) << ", " << i->r(1) << ", "<< i->r(2) << ", F: "<< i->F(0) << ", "<< i->F(1) << ", "<< i->F(2) << std::endl;
		}
		if (containerTest->requiresForceExchange()) {
			domainDecompositionTest->exchangeForces(containerTest, _domain);
		}
//		const ParticleIterator& begin2 = containerTest->iteratorBegin();
//		const ParticleIterator& end2 = containerTest->iteratorEnd();
//		for (ParticleIterator i = begin2; i != end2; ++i) {
//			std::cout << "r: " << i->r(0) << ", " << i->r(1) << ", " << i->r(2) << ", F: " << i->F(0) << ", " << i->F(1)
//					<< ", " << i->F(2) << std::endl;
//		}
	}
	//------------------------------------------------------------
	container->deleteOuterParticles();
	containerTest->deleteOuterParticles();

	// Compare calculated forces
	{
		const ParticleIterator begin = container->iteratorBegin();
		const ParticleIterator end = container->iteratorEnd();
		const ParticleIterator beginHS = containerTest->iteratorBegin();
		const ParticleIterator endHS = containerTest->iteratorEnd();
		auto j = beginHS;
		for (auto i = begin; i != end; ++i, ++j) {
			CPPUNIT_ASSERT_EQUAL(j->id(), i->id());
#ifdef MARDYN_SPSP
			double delta=1e-6;
#else
			double delta=1e-9;
#endif
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Forces differ", i->F(0), j->F(0), fabs(delta * i->F(0)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Forces differ", i->F(1), j->F(1), fabs(delta * i->F(1)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Forces differ", i->F(2), j->F(2), fabs(delta * i->F(2)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Virials differ", i->Vi(0), j->Vi(0), fabs(delta * j->Vi(0)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Virials differ", i->Vi(1), j->Vi(1), fabs(delta * j->Vi(1)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Virials differ", i->Vi(2), j->Vi(2), fabs(delta * j->Vi(2)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Rotational moments differ", i->M(0), j->M(0), fabs(delta * i->M(0)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Rotational moments differ", i->M(1), j->M(1), fabs(delta * i->M(1)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Rotational moments differ", i->M(2), j->M(2), fabs(delta * i->M(2)));
		}
	}
	//------------------------------------------------------------
	// Cleanup
	//------------------------------------------------------------

	delete domainDecompositionFS;
	delete domainDecompositionTest;
	delete cellProc;
	delete cellProc2;
	delete container;
	delete containerTest;

}
