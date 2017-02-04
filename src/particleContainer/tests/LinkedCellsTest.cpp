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

TEST_SUITE_REGISTRATION(LinkedCellsTest);

LinkedCellsTest::LinkedCellsTest() {

}

LinkedCellsTest::~LinkedCellsTest() {
}

void LinkedCellsTest::testUpdateAndDeleteOuterParticlesFilename(const char * filename, double cutoff) {
	// original pointer will be deleted by tearDown()
	_domainDecomposition = new DomainDecompBase();

	LinkedCells* container = static_cast<LinkedCells*> (initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff));
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
	Molecule dummyMolecule1(1, &_components[0], 1.0,1.0,1.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule2(2, &_components[0], 2.0,2.0,2.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule3(3, &_components[0], 3.0,3.0,3.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule4(4, &_components[0], 5.1,5.1,5.1,0,0,0, 0, 0, 0, 0, 0, 0, 0);

	LinkedCells LC;
	MoleculeIterator molIt;

	std::vector<ParticleCell> & cells = LC._cells;

	cells.resize(5);

	// cell[0] is empty

	// cell[1] contains 1,2,3
	{
		double l[3] = {0., 0., 0.}, u[3] = {5., 5., 5.};
		cells[1].setBoxMin(l);
		cells[1].setBoxMax(u);
	}
	cells[1].addParticle(dummyMolecule1);
	cells[1].addParticle(dummyMolecule2);
	cells[1].addParticle(dummyMolecule3);

	// cell[2] is empty

	// cell[3] contains 4
	{
		double l[3] = {5., 5., 5.}, u[3] = {5.5, 5.5, 5.5};
		cells[3].setBoxMin(l);
		cells[3].setBoxMax(u);
	}
	cells[3].addParticle(dummyMolecule4);

	// cell[4] is empty

	// NEXTNONEMPTY: start at particle 3, and arrive at particle 4

	LC._cellIterator = cells.begin() + 1;
	molIt = LC.nextNonEmptyCell();
	ASSERT_TRUE_MSG("nextNonEmpty::_cellIterator", LC._cellIterator == LC._cells.begin() + 3);
	ASSERT_TRUE_MSG("nextNonEmpty::_particleIterator", LC.current()->id() == 4ul);
	ASSERT_TRUE_MSG("nextNonEmpty::return", molIt->id() == 4ul);

	// BEGIN:
	molIt = LC.begin();
	ASSERT_TRUE_MSG("begin()", molIt->id() == 1ul);
	ASSERT_TRUE_MSG("end()", molIt != LC.end());

	// NEXT:
	molIt = LC.next();
	ASSERT_TRUE_MSG("next() within cell", molIt->id() == 2ul);
	ASSERT_TRUE_MSG("end()", molIt != LC.end());
	molIt = LC.next();
	ASSERT_TRUE_MSG("next() within cell", molIt->id() == 3ul);
	ASSERT_TRUE_MSG("end()", molIt != LC.end());
	molIt = LC.next();
	ASSERT_TRUE_MSG("next() across cells", molIt->id() == 4ul);
	ASSERT_TRUE_MSG("end()", molIt != LC.end());
	molIt = LC.next();
	ASSERT_TRUE_MSG("next() arrive at end()", molIt == LC.end());

	// DELETECURRENT:
	molIt = LC.begin();

	molIt = LC.deleteCurrent();
	ASSERT_EQUAL_MSG("delete() within cell", 3ul, molIt->id()); // 3 copied in place of 1
	molIt = LC.deleteCurrent();
	ASSERT_TRUE_MSG("delete() within cell", molIt->id() == 2ul); // 2 copied in place of 3
	molIt = LC.deleteCurrent();
	ASSERT_TRUE_MSG("delete() across cells", molIt->id() == 4ul); // cell 1 became empty, we advanced to cell 3
	molIt = LC.deleteCurrent();
	ASSERT_TRUE_MSG("delete() last", molIt == LC.end()); // cell 4 became empty, we arrived at end()
}

void LinkedCellsTest::testParticleIteratorBeginNextEndParticleIteratorSequential() {
	// NOTE: we do not open an OpenMP parallel region!
	// Hence, this test is always executed sequentially!

	Molecule dummyMolecule1(1, &_components[0], 1.0,1.0,1.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule2(2, &_components[0], 2.0,2.0,2.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule3(3, &_components[0], 3.0,3.0,3.0,0,0,0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule4(4, &_components[0], 5.1,5.1,5.1,0,0,0, 0, 0, 0, 0, 0, 0, 0);

	LinkedCells LC;

	std::vector<ParticleCell> & cells = LC._cells;

	cells.resize(5);

	// cell[0] is empty

	// cell[1] contains 1,2,3
	{
		double l[3] = {0., 0., 0.}, u[3] = {5., 5., 5.};
		cells[1].setBoxMin(l);
		cells[1].setBoxMax(u);
	}
	cells[1].addParticle(dummyMolecule1);
	cells[1].addParticle(dummyMolecule2);
	cells[1].addParticle(dummyMolecule3);

	// cell[2] is empty

	// cell[3] contains 4
	{
		double l[3] = {5., 5., 5.}, u[3] = {5.5, 5.5, 5.5};
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
			CellProcessor(0., 0.), _num_cells(num_cells) {
		_cellProcessCount.assign(num_cells, 0);
		_cellPairProcessCount.assign(num_cells * num_cells, 0);
	}
	virtual void initTraversal() {
	}

	virtual void preprocessCell(ParticleCell& /*cell*/) {
	}

	virtual void processCellPair(ParticleCell& cell1, ParticleCell& cell2) {
		_cellPairProcessCount[_num_cells * cell1.getCellIndex() + cell2.getCellIndex()] += sign;
		_cellPairProcessCount[_num_cells * cell2.getCellIndex() + cell1.getCellIndex()] += sign;  // newton 3
	}

	virtual void processCell(ParticleCell& cell) {
		_cellProcessCount[cell.getCellIndex()] += sign;
	}

	virtual double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) {
		return 0.;
	}

	virtual int countNeighbours(Molecule* /*m1*/, ParticleCell& /*cell2*/, double /*RR*/) {
		return 0;
	}

	virtual void postprocessCell(ParticleCell& /*cell*/) {
	}

	virtual void endTraversal() {

	}

	void checkZero() {
		for (int i : _cellProcessCount) {
			ASSERT_EQUAL(i, 0);
		}
		for (int i : _cellPairProcessCount) {
			ASSERT_EQUAL(i, 0);
		}
	}

	void checkOnlyInner(LinkedCells* container) {
		for (int i = 0; i < _num_cells; i++) {
			if (_cellProcessCount[i] != 0) {
				ASSERT_TRUE(container->getCell(i).isInnerCell());
			}
		}
		for (int i = 0; i < _num_cells; i++) {
			for (int j = 0; j < _num_cells; j++) {
				if (_cellPairProcessCount[i * _num_cells + j] != 0) {
					ASSERT_TRUE(container->getCell(i).isInnerCell() && container->getCell(j).isInnerCell());
				}
			}
		}
	}

	void inverseSign() {
		sign *= -1;
	}
private:
	std::vector<int> _cellProcessCount;
	std::vector<int> _cellPairProcessCount;
	int _num_cells;
	int sign = 1;
};

void LinkedCellsTest::testTraversalMethods() {
	const char* filename = "VectorizationMultiComponentMultiPotentials.inp";
	ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 25.);
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
