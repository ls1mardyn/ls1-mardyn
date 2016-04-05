/*
 * LinkedCellsTest.cpp
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#include "LinkedCellsTest.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

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
	cells[1].addParticle(new Molecule(dummyMolecule1));
	cells[1].addParticle(new Molecule(dummyMolecule2));
	cells[1].addParticle(new Molecule(dummyMolecule3));

	// cell[2] is empty

	// cell[3] contains 4
	cells[3].addParticle(new Molecule(dummyMolecule4));

	// cell[4] is empty

	// NEXTNONEMPTY: start at particle 3, and arrive at particle 4

	LC._cellIterator = cells.begin() + 1;
	molIt = LC.nextNonEmptyCell();
	ASSERT_TRUE_MSG("nextNonEmpty::_cellIterator", LC._cellIterator == LC._cells.begin() + 3);
	ASSERT_TRUE_MSG("nextNonEmpty::_particleIterator", (*LC._particleIterator)->id() == 4ul);
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
	ASSERT_TRUE_MSG("delete() within cell", molIt->id() == 3ul); // 3 copied in place of 1
	molIt = LC.deleteCurrent();
	ASSERT_TRUE_MSG("delete() within cell", molIt->id() == 2ul); // 2 copied in place of 3
	molIt = LC.deleteCurrent();
	ASSERT_TRUE_MSG("delete() across cells", molIt->id() == 4ul); // cell 1 became empty, we advanced to cell 3
	molIt = LC.deleteCurrent();
	ASSERT_TRUE_MSG("delete() last", molIt == LC.end()); // cell 4 became empty, we arrived at end()
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
