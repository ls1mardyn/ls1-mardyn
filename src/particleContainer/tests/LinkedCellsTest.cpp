/*
 * LinkedCellsTest.cpp
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#include "LinkedCellsTest.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/adapter/CellProcessor.h"
#include <vector>

#include "ensemble/PressureGradient.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/LinkedCellTraversals/HalfShellTraversal.h"

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
}

void LinkedCellsTest::testUpdateAndDeleteOuterParticlesH2O() {
	testUpdateAndDeleteOuterParticlesFilename("H20_NaBr_0.01_T_293.15.inp", 27.0);
}

void LinkedCellsTest::testUpdateAndDeleteOuterParticles8Particles() {
	testUpdateAndDeleteOuterParticlesFilename("LinkedCells.inp", 1.0);
}

void LinkedCellsTest::testMoleculeBeginNextEndDeleteCurrent() {
	// NOTE: we do not open an OpenMP parallel region!
	// Hence, this test is always executed sequentially!

	Molecule dummyMolecule1(1, &_components[0], 1.0, 1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule2(2, &_components[0], 2.0, 2.0, 2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule3(3, &_components[0], 3.0, 3.0, 3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule4(4, &_components[0], 5.1, 5.1, 5.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	LinkedCells LC;
	ParticleIterator molIt;

	std::vector<ParticleCell> & cells = LC._cells;

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
	// NOTE: we do not open an OpenMP parallel region!
	// Hence, this test is always executed sequentially!

	Molecule dummyMolecule1(1, &_components[0], 1.0, 1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule2(2, &_components[0], 2.0, 2.0, 2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule3(3, &_components[0], 3.0, 3.0, 3.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule4(4, &_components[0], 5.1, 5.1, 5.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

	LinkedCells LC;

	std::vector<ParticleCell> & cells = LC._cells;

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
		CellProcessor(0., 0.), _cellProcessCount(num_cells, 0), _cellPairProcessCount(num_cells),
		_num_cells(num_cells) {
	}
	virtual void initTraversal() {
	}

	virtual void preprocessCell(ParticleCell& /*cell*/) {
	}

	virtual void processCellPair(ParticleCell& cell1, ParticleCell& cell2) {
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
					ASSERT_TRUE(container->getCellReference(i).isInnerCell() && container->getCellReference(pair.first).isInnerCell());
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

}

void LinkedCellsTest::testHalfShell() {
	//TODO: ___Extract to separate test class
	//------------------------------------------------------------
	// Setup
	//------------------------------------------------------------

	auto domainDecomposition = new DomainDecompBase();
	auto cutoff = 1;
	auto filename = "LinkedCellsHS.inp";

	LinkedCells* containerHS = dynamic_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
			filename, cutoff));
	containerHS->_traversalSelected = LinkedCells::Traversal::HS;
	containerHS->_traversal = nullptr;
	containerHS->initializeTraversal();

	LinkedCells* container = dynamic_cast<LinkedCells*>(initializeFromFile(ParticleContainerFactory::LinkedCell,
			filename, cutoff));

	auto vectorizedCellProcessor = new VectorizedCellProcessor(*_domain, cutoff, cutoff);

	//------------------------------------------------------------
	// Prepare molecule containers
	//------------------------------------------------------------

	container->update();
	containerHS->update();

	bool forceRebalancing = false;
	domainDecomposition->balanceAndExchange(forceRebalancing, container, _domain);
	domainDecomposition->balanceAndExchange(forceRebalancing, containerHS, _domain);

	container->updateMoleculeCaches();
	containerHS->updateMoleculeCaches();

	//------------------------------------------------------------
	// Do calculation with HS
	//------------------------------------------------------------
	{
		containerHS->traverseCells(*vectorizedCellProcessor);

		// calculate forces
		const ParticleIterator& begin = containerHS->iteratorBegin();
		const ParticleIterator& end = containerHS->iteratorEnd();
		for (ParticleIterator i = begin; i != end; ++i) {
			i->calcFM();
		}

		domainDecomposition->exchangeForces(containerHS, _domain);
	}
	//------------------------------------------------------------
	// Do calculation with FS
	//------------------------------------------------------------
	{
		container->traverseCells(*vectorizedCellProcessor);

		// calculate forces
		const ParticleIterator begin = container->iteratorBegin();
		const ParticleIterator end = container->iteratorEnd();
		for (auto i = begin; i != end; ++i) {
			i->calcFM();
		}
	}
	//------------------------------------------------------------
	container->deleteOuterParticles();
	containerHS->deleteOuterParticles();
	// Compare calculated forces
	{
		const ParticleIterator begin = container->iteratorBegin();
		const ParticleIterator end = container->iteratorEnd();
		const ParticleIterator beginHS = containerHS->iteratorBegin();
		const ParticleIterator endHS = containerHS->iteratorEnd();
		auto j = beginHS;
		for (auto i = begin; i != end; ++i, ++j) {

			CPPUNIT_ASSERT_EQUAL(j->id(), i->id());
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Forces differ", i->F(0), j->F(0), fabs(1e-7*i->F(0)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Forces differ", i->F(1), j->F(1), fabs(1e-7*i->F(1)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Forces differ", i->F(2), j->F(2), fabs(1e-7*i->F(2)));

			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Virials differ", i->Vi(0), j->Vi(0), fabs(1e-7*j->Vi(0)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Virials differ", i->Vi(1), j->Vi(1), fabs(1e-7*j->Vi(1)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Virials differ", i->Vi(2), j->Vi(2), fabs(1e-7*j->Vi(2)));

			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Rotational moments differ", i->M(0), j->M(0), fabs(1e-7*i->M(0)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Rotational moments differ", i->M(1), j->M(1), fabs(1e-7*i->M(1)));
			CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE("Rotational moments differ", i->M(2), j->M(2), fabs(1e-7*i->M(2)));
		}
	}

	//------------------------------------------------------------
	// Cleanup
	//------------------------------------------------------------

	delete domainDecomposition;
	delete vectorizedCellProcessor;
}
