/**
 * @file ODFCellProcessor.cpp
 * @author F. Gratl
 * @date 10/11/19
 */

#include "ODFCellProcessor.h"
#include "particleContainer/ParticleCell.h"
#include "io/ODF.h"

ODFCellProcessor::ODFCellProcessor(const double cutoffRadius, ODF *const odf, const std::array<double, 3> &simBoxSize)
	: CellProcessor(cutoffRadius, cutoffRadius), _odf(odf), _simBoxSize(simBoxSize) {}

void ODFCellProcessor::initTraversal() {}
void ODFCellProcessor::preprocessCell(ParticleCell &cell) {}
void ODFCellProcessor::processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) {
	auto begin1 = cell1.iterator();
	auto begin2 = cell2.iterator();

	if (sumAll or not cell2.isHaloCell() or cell1.isInnerCell() or cell1.getCellIndex() < cell2.getCellIndex()) {
		// loop over all particles in the cell
		for (auto it1 = begin1; it1.isValid(); ++it1) {
			Molecule &molecule1 = *it1;

			if (molecule1.numDipoles() != 1) {
				continue;
			}

			auto orientationVector1 = calcOrientationVector(molecule1);

			for (auto it2 = begin2; it2.isValid(); ++it2) {
				Molecule &molecule2 = *it2;

				if (molecule2.numDipoles() != 1) {
					continue;
				}

				double dummy[3];
				double dd = molecule2.dist2(molecule1, dummy);
				if (dd < _cutoffRadiusSquare) {
					_odf->calculateOrientation(_simBoxSize, molecule1, molecule2, orientationVector1);
				}
			}
		}
	}
}
void ODFCellProcessor::processCell(ParticleCell &cell) {
	if (cell.isInnerCell() || cell.isBoundaryCell()) {
		auto begin = cell.iterator();

		for (auto it1 = begin; it1.isValid(); ++it1) {
			Molecule &molecule1 = *it1;

			if (molecule1.numDipoles() != 1) {
				continue;
			}

			auto orientationVector1 = calcOrientationVector(molecule1);

			// get second molecule
			auto it2 = it1;
			++it2;
			for (; it2.isValid(); ++it2) {
				Molecule &molecule2 = *it2;

				mardyn_assert(&molecule1 != &molecule2);

				if (molecule2.numDipoles() != 1) {
					continue;
				}

				double dummy[3];
				double dd = molecule2.dist2(molecule1, dummy);

				if (dd < _cutoffRadiusSquare) {
					_odf->calculateOrientation(_simBoxSize, molecule1, molecule2, orientationVector1);
				}
			}
		}
	}
}

double ODFCellProcessor::processSingleMolecule(Molecule *m1, ParticleCell &cell2) { return 0; }
void ODFCellProcessor::postprocessCell(ParticleCell &cell) {}
void ODFCellProcessor::endTraversal() {}

std::array<double, 3> ODFCellProcessor::calcOrientationVector(const Molecule &molecule) {
	std::array<double, 4> q1{molecule.q().qw(), molecule.q().qx(), molecule.q().qy(), molecule.q().qz()};
	std::array<double, 3> orientationVector = {2 * (q1[1] * q1[3] + q1[0] * q1[2]), 2 * (q1[2] * q1[3] - q1[0] * q1[1]),
											   1 - 2 * (q1[1] * q1[1] + q1[2] * q1[2])};
	return orientationVector;
}
