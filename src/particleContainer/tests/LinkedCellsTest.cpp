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
