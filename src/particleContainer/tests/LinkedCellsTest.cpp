/*
 * LinkedCellsTest.cpp
 *
 * @Date: 03.05.2011
 * @Author: eckhardw
 */

#include "LinkedCellsTest.h"
#include "parallel/DomainDecompBase.h"

TEST_SUITE_REGISTRATION(LinkedCellsTest);

LinkedCellsTest::LinkedCellsTest() {

}

LinkedCellsTest::~LinkedCellsTest() {
}

void LinkedCellsTest::testUpdateAndDeleteOuterParticlesH2O() {

	// original pointer will be deleted by tearDown()
	_domainDecomposition = new DomainDecompBase();

	LinkedCells* container = static_cast<LinkedCells*> (initializeFromFile(ParticleContainerFactory::LinkedCell, "H20_NaBr_0.01_T_293.15.inp", 27.0));
	int numMols = container->getNumberOfParticles();

	_domainDecomposition->exchangeMolecules(container, _domain);
	container->deleteOuterParticles();

	int newNumMols = container->getNumberOfParticles();
	ASSERT_EQUAL(numMols, newNumMols);

	delete _domainDecomposition;
}
