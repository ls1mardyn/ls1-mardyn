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

#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"
#include "particleContainer/LinkedCellTraversals/HalfShellTraversal.h"
#include "particleContainer/TraversalTuner.h"

#include "utils/arrayMath.h"

TEST_SUITE_REGISTRATION(LinkedCellsTest);

LinkedCellsTest::LinkedCellsTest() {

}

LinkedCellsTest::~LinkedCellsTest() {
}

void LinkedCellsTest::testRegionIterator() {
	// can't be const bc of LinkedCells constructor...
	std::array<double, 3> boxMin = {0.0, 0.0, 0.0};
	std::array<double, 3> boxMax = {10.0, 10.0, 10.0};
	std::array<double, 3> boxLength = arrayMath::sub(boxMax, boxMin);
	const double cutoff = 1.;
	LinkedCells linkedCells(boxMin.data(), boxMax.data(), cutoff);

	size_t numMolecules = 0;
	// place one molecule in the center of every cell
	for (double z = boxMin[2]; z < boxMax[2]; ++z) {
		for (double y = boxMin[1]; y < boxMax[1]; ++y) {
			for (double x = boxMin[0]; x < boxMax[0]; ++x) {
				Molecule m(numMolecules++, &_components[0], x + cutoff / 2, y + cutoff / 2, z + cutoff / 2, 0., 0., 0., 0., 0.,
						   0., 0., 0., 0., 0.);
				linkedCells.addParticle(m);
			}
		}
	}

	// Test over full domain
	{
		size_t particlesFound = 0;
		const std::array<double, 3> iterRegionStart = boxMin;
		const std::array<double, 3> iterRegionEnd = boxMax;
		for (auto regionIter = linkedCells.regionIterator(iterRegionStart.data(), iterRegionEnd.data(), ParticleIterator::ALL_CELLS);
			 regionIter.isValid(); ++regionIter) {
			++particlesFound;
		}
		ASSERT_EQUAL_MSG("RegionIterator over the whole domain did not find all particles.",
				numMolecules, particlesFound);
	}
	// Test outside the domain
	{
		size_t particlesFound = 0;
		const std::array<double, 3> iterRegionStart = {boxMax[0] + 10., boxMin[1] - 1., boxMin[2] - 1.};
		const std::array<double, 3> iterRegionEnd = {boxMax[0] + 20., boxMax[1] + 1., boxMax[2] + 1.};
		for (auto regionIter = linkedCells.regionIterator(iterRegionStart.data(), iterRegionEnd.data(), ParticleIterator::ALL_CELLS);
			 regionIter.isValid(); ++regionIter) {
			++particlesFound;
		}
		ASSERT_EQUAL_MSG("RegionIterator outside the domain found particles, but it shouldn't.",
				0ul, particlesFound);
	}
	// Test partially outside the domain
	{
		size_t particlesFound = 0;
		const std::array<double, 3> iterRegionStart = {boxMin[0] + 0.5*boxLength[0], boxMin[1] - 1., boxMin[2] - 1.};
		const std::array<double, 3> iterRegionEnd = {boxMax[0] + 20., boxMax[1] + 1., boxMax[2] + 1.};
		for (auto regionIter = linkedCells.regionIterator(iterRegionStart.data(), iterRegionEnd.data(), ParticleIterator::ALL_CELLS);
			 regionIter.isValid(); ++regionIter) {
			++particlesFound;
		}
		ASSERT_EQUAL_MSG("RegionIterator outside the domain found particles, but it shouldn't.",
				numMolecules/2, particlesFound);
	}
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
	// NOTE: we do not open an OpenMP parallel region!
	// Hence, this test is always executed sequentially!

	Molecule dummyMolecule1(1, &_components[0], 0.1, 0.1, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule2(2, &_components[0], 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule3(3, &_components[0], 0.3, 0.3, 0.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	Molecule dummyMolecule4(4, &_components[0], 0.1, 1.5, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);


	double bBoxMin[3] = {0.0, 0.0, 0.0};
	double bBoxMax[3] = {2.0, 2.0, 2.0};
	double cutoffRadius = 1.0;
	LinkedCells LC(bBoxMin, bBoxMax, cutoffRadius);

	// some empty cells

	// particles 1,2,3 in same cell
	LC.addParticle(dummyMolecule1);
	LC.addParticle(dummyMolecule2);
	LC.addParticle(dummyMolecule3);

	// some empty cells

	// particle 4
	LC.addParticle(dummyMolecule4);

	// some empty cells

	// BEGIN:
	auto molIt = LC.iterator(ParticleIterator::ALL_CELLS);
	ASSERT_TRUE_MSG("begin()", molIt->getID() == 1ul);
	ASSERT_TRUE_MSG("end()", molIt.isValid());

	// NEXT:
	++molIt;
	ASSERT_TRUE_MSG("next() within cell", molIt->getID() == 2ul);
	ASSERT_TRUE_MSG("end()", molIt.isValid());
	++molIt;
	ASSERT_TRUE_MSG("next() within cell", molIt->getID() == 3ul);
	ASSERT_TRUE_MSG("end()", molIt.isValid());
	++molIt;
	ASSERT_TRUE_MSG("next() across cells", molIt->getID() == 4ul);
	ASSERT_TRUE_MSG("end()", molIt.isValid());
	++molIt;
	ASSERT_TRUE_MSG("next() arrive at end()", not molIt.isValid());

	// DELETECURRENT:
	molIt = LC.iterator(ParticleIterator::ALL_CELLS);

	molIt.deleteCurrentParticle();
	++molIt;
	ASSERT_EQUAL_MSG("delete() within cell", 3ul, molIt->getID()); // 3 copied in place of 1
	molIt.deleteCurrentParticle();
	++molIt;
	ASSERT_TRUE_MSG("delete() within cell", molIt->getID() == 2ul); // 2 copied in place of 3
	molIt.deleteCurrentParticle();
	++molIt;
	ASSERT_TRUE_MSG("delete() across cells", molIt->getID() == 4ul); // cell 1 became empty, we advanced to cell 3
	molIt.deleteCurrentParticle();
	++molIt;
	ASSERT_TRUE_MSG("delete() last", not molIt.isValid()); // cell 4 became empty, we arrived at end()
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
		_cellPairProcessCount[cell2.getCellIndex()][cell1.getCellIndex()] += sign;	// newton 3
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
	int* boxWidthInNumCells = dynamic_cast<LinkedCells*>(container)->getBoxWidthInNumCells();
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
//	//	Calculate forces for FS and HS and compare
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
//	//	Calculate forces for FS and MP and compare
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

void LinkedCellsTest::testFullShellMPIDirectPP() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "indirect", "hs");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::C08, 1, "direct-pp", "fs");
}

void LinkedCellsTest::testFullShellMPIDirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "indirect", "hs");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::C08, 1, "direct", "fs");
}

void LinkedCellsTest::testHalfShellMPIIndirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "indirect", "hs");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "indirect", "hs");
}

void LinkedCellsTest::testHalfShellMPIDirectPP() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "direct", "hs");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "direct-pp", "hs");
}

void LinkedCellsTest::testHalfShellMPIDirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "direct", "hs");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::HS, 1, "direct", "hs");
}

void LinkedCellsTest::testMidpointMPIIndirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "indirect", "mp");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "indirect", "mp");
}

void LinkedCellsTest::testMidpointMPIDirectPP() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "direct", "mp");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "direct-pp", "mp");
}

void LinkedCellsTest::testMidpointMPIDirect() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "direct", "mp");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "direct", "mp");
}

void LinkedCellsTest::testEighthShellMPIDirectPP() {
//	doForceComparisonTest("simple-lj.inp", TraversalTuner < ParticleCell > ::traversalNames::MP, 2, "direct", "mp");
	doForceComparisonTest("simple-lj-tiny.inp", TraversalTuner < ParticleCell > ::traversalNames::C08ES, 2, "direct-pp", "es");
}

void LinkedCellsTest::testCellBorderAndFlagManager() {
	long int cellIndex;
	double cellBoxMin[3], cellBoxMax[3];

	double bMin[3] = {0.1, 0.2, 0.3};
	double bMax[3] = {5.1, 6.1, 7.3};
	double cutoff = 0.7;
	LinkedCells LC(bMin, bMax, cutoff);

	// this contains all of the necessary assert-statements
	for (int iz = 0; iz < LC._cellsPerDimension[2]; ++iz) {
		cellBoxMin[2] = iz * LC._cellLength[2] + LC._haloBoundingBoxMin[2];
		cellBoxMax[2] = (iz + 1) * LC._cellLength[2] + LC._haloBoundingBoxMin[2];
		if (iz == 0) {	// make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
			cellBoxMax[2] = LC._boundingBoxMin[2];
		} else if (iz == 1) {// make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
			cellBoxMin[2] = LC._boundingBoxMin[2];
		} else if (iz == LC._cellsPerDimension[2] - 2) { // make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
			cellBoxMax[2] = LC._boundingBoxMax[2];
		} else if (iz == LC._cellsPerDimension[2] - 1) { // make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
			cellBoxMin[2] = LC._boundingBoxMax[2];
			cellBoxMax[2] = LC._haloBoundingBoxMax[2];
		}
		for (int iy = 0; iy < LC._cellsPerDimension[1]; ++iy) {
			cellBoxMin[1] = iy * LC._cellLength[1] + LC._haloBoundingBoxMin[1];
			cellBoxMax[1] = (iy + 1) * LC._cellLength[1] + LC._haloBoundingBoxMin[1];
			if (iy == 0) { // make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
				cellBoxMax[1] = LC._boundingBoxMin[1];
			} else if (iy == 1) {// make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
				cellBoxMin[1] = LC._boundingBoxMin[1];
			} else if (iy == LC._cellsPerDimension[1] - 2) { // make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
				cellBoxMax[1] = LC._boundingBoxMax[1];
			} else if (iy == LC._cellsPerDimension[1] - 1) { // make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
				cellBoxMin[1] = LC._boundingBoxMax[1];
				cellBoxMax[1] = LC._haloBoundingBoxMax[1];
			}
			for (int ix = 0; ix < LC._cellsPerDimension[0]; ++ix) {
				cellBoxMin[0] = ix * LC._cellLength[0] + LC._haloBoundingBoxMin[0];
				cellBoxMax[0] = (ix + 1) * LC._cellLength[0] + LC._haloBoundingBoxMin[0];
				if (ix == 0) { // make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
					cellBoxMax[0] = LC._boundingBoxMin[0];
				} else if (ix == 1) {// make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
					cellBoxMin[0] = LC._boundingBoxMin[0];
				} else if (ix == LC._cellsPerDimension[0] - 2) { // make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
					cellBoxMax[0] = LC._boundingBoxMax[0];
				} else if (ix == LC._cellsPerDimension[0] - 1) { // make sure, that the cells span the whole domain... for iz=0 this is already implicitly done
					cellBoxMin[0] = LC._boundingBoxMax[0];
					cellBoxMax[0] = LC._haloBoundingBoxMax[0];
				}
				cellIndex = LC.cellIndexOf3DIndex(ix, iy, iz);
				ParticleCell & cell = LC._cells[cellIndex];
				cell.setCellIndex(cellIndex); //set the index of the cell to the index of it...

				cell.setBoxMin(cellBoxMin);
				cell.setBoxMax(cellBoxMax);
				if (ix < LC._haloWidthInNumCells[0] ||
					iy < LC._haloWidthInNumCells[1] ||
					iz < LC._haloWidthInNumCells[2] ||
					ix >= LC._cellsPerDimension[0] - LC._haloWidthInNumCells[0] ||
					iy >= LC._cellsPerDimension[1] - LC._haloWidthInNumCells[1] ||
					iz >= LC._cellsPerDimension[2] - LC._haloWidthInNumCells[2]) {

					cell.assignCellToHaloRegion();
				}
				else{
					if (ix < 2 * LC._haloWidthInNumCells[0] ||
						iy < 2 * LC._haloWidthInNumCells[1] ||
						iz < 2 * LC._haloWidthInNumCells[2] ||
						ix >= LC._cellsPerDimension[0] - 2 * LC._haloWidthInNumCells[0] ||
						iy >= LC._cellsPerDimension[1] - 2 * LC._haloWidthInNumCells[1] ||
						iz >= LC._cellsPerDimension[2] - 2 * LC._haloWidthInNumCells[2]) {

						cell.assignCellToBoundaryRegion();
					}
					else{
						if (ix < 3 * LC._haloWidthInNumCells[0] ||
							iy < 3 * LC._haloWidthInNumCells[1] ||
							iz < 3 * LC._haloWidthInNumCells[2] ||
							ix >= LC._cellsPerDimension[0] - 3 * LC._haloWidthInNumCells[0] ||
							iy >= LC._cellsPerDimension[1] - 3 * LC._haloWidthInNumCells[1] ||
							iz >= LC._cellsPerDimension[2] - 3 * LC._haloWidthInNumCells[2]) {

							cell.assignCellToInnerRegion();
						}
						else {
							cell.assignCellToInnerMostAndInnerRegion();
						}
					}
				}
			}
		}
	}

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
	domainDecompositionFS->setCommunicationScheme("indirect", "fs");
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
	domainDecompositionFS->initCommunicationPartners(cutoff, _domain, container);
	domainDecompositionTest->initCommunicationPartners(cutoff, _domain, containerTest);
#endif

	//------------------------------------------------------------
	//	Calculate forces for FS and TestTraversal and compare
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
	//std::cout << std::endl<<"reference:"<<std::endl;
	{
		container->traverseCells(*cellProc);
		// calculate forces
		const auto begin = container->iterator(ParticleIterator::ALL_CELLS);
		for (auto i = begin; i.isValid(); ++i) {
			i->calcFM();
	//		std::cout << "r: " << i->r(0) << ", " << i->r(1) << ", "<< i->r(2) << ", F: "<< i->F(0) << ", "<< i->F(1) << ", "<< i->F(2) << std::endl;
		}
	}
	//std::cout << std::endl <<"test:"<< std::endl;
	//------------------------------------------------------------
	// Do calculation with TestTraversal
	//------------------------------------------------------------
	{
		containerTest->traverseCells(*cellProc2);
		// calculate forces

	//	std::cout << "pre force exchange:" << std::endl;
		for (auto i = containerTest->iterator(ParticleIterator::ALL_CELLS); i.isValid(); ++i) {
			i->calcFM();
	//		std::cout << "r: " << i->r(0) << ", " << i->r(1) << ", "<< i->r(2) << ", F: "<< i->F(0) << ", "<< i->F(1) << ", "<< i->F(2) << std::endl;
		}
		if (containerTest->requiresForceExchange()) {
			domainDecompositionTest->exchangeForces(containerTest, _domain);
		}
		std::cout << "after force exchange:"<< std::endl;

		for (auto i = containerTest->iterator(ParticleIterator::ALL_CELLS); i.isValid(); ++i) {
	//		std::cout << "r: " << i->r(0) << ", " << i->r(1) << ", " << i->r(2) << ", F: " << i->F(0) << ", " << i->F(1)
	//				<< ", " << i->F(2) << std::endl;
		}
	}
	//------------------------------------------------------------
	container->deleteOuterParticles();
	containerTest->deleteOuterParticles();

	// Compare calculated forces
	{
		const auto begin = container->iterator(ParticleIterator::ALL_CELLS);
		const auto beginHS = containerTest->iterator(ParticleIterator::ALL_CELLS);
		auto j = beginHS;
		for (auto i = begin; i.isValid(); ++i, ++j) {
			CPPUNIT_ASSERT_EQUAL(j->getID(), i->getID());
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
