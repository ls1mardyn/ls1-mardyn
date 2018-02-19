/*
 * RDFCellProcessor.cpp
 *
 * @Date: 19.02.18
 * @Author: tchipevn
 */

#include "RDFCellProcessor.h"
#include "particleContainer/ParticleCell.h"
#include "io/RDF.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include <vector>

using namespace std;
using namespace Log;

void RDFCellProcessor::processCell(ParticleCell& cell) {
	if (cell.isInnerCell() || cell.isBoundaryCell()) {
		SingleCellIterator begin = cell.iteratorBegin();
		SingleCellIterator end = cell.iteratorEnd();

		for (SingleCellIterator it1 = begin; it1 != end; ++it1) {
			Molecule& molecule1 = *it1;

			SingleCellIterator it2 = it1;
			++it2;
			for (; it2 != end; ++it2) {
				Molecule& molecule2 = *it2;

				mardyn_assert(&molecule1 != &molecule2);
				double dummy[3];
				double dd = molecule2.dist2(molecule1, dummy);

				if (dd < _cutoffRadiusSquare) {
					_rdf->observeRDF(molecule1, molecule2, dd);
				}
			}
		}
	} // if (isInnerCell())
}

void RDFCellProcessor::processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll /* = false */) {
	SingleCellIterator begin1 = cell1.iteratorBegin();
	SingleCellIterator end1 = cell1.iteratorEnd();
	SingleCellIterator begin2 = cell2.iteratorBegin();
	SingleCellIterator end2 = cell2.iteratorEnd();

	if(sumAll) { // sumAll - moleculesAt is gone, use SingleCellIterator now ?

		// loop over all particles in the cell
		for (SingleCellIterator it1 = begin1; it1 != end1; ++it1) {
			Molecule& molecule1 = *it1; 
			for (SingleCellIterator it2 = begin2; it2 != end2; ++it2) {
				Molecule& molecule2 = *it2; 

				double dummy[3];
				double dd = molecule2.dist2(molecule1, dummy);
				if (dd < _cutoffRadiusSquare) {
                    _rdf->observeRDF(molecule1, molecule2, dd);
				}
			}

		}
	} else { // sumHalf
		if (cell1.isInnerCell()) {//no cell is halo
			// loop over all particles in the cell


			for (SingleCellIterator it1 = begin1; it1 != end1; ++it1) {
				Molecule& molecule1 = *it1;

				for (SingleCellIterator it2 = begin2; it2 != end2; ++it2) {
					Molecule& molecule2 = *it2;

					double dummy[3];
					double dd = molecule2.dist2(molecule1, dummy);
					if (dd < _cutoffRadiusSquare) {
	                    _rdf->observeRDF(molecule1, molecule2, dd);
					}
				}

			}
		} // inner cell

		if (cell1.isBoundaryCell()) {//first cell is boundary
			//if (cell2.isHaloCell() && ! molecule1.isLessThan(molecule2)) {//boundary <-> halo: using macroscopic boundary condition
			if (cell2.isHaloCell() && ! (cell1.getCellIndex()<cell2.getCellIndex())) {//boundary <-> halo: not using macroscopic boundary condition, instead cell indices are compared.
				/* Do not sum up values twice. */
				return;
			}

			for (SingleCellIterator it1 = begin1; it1 != end1; ++it1) {
				Molecule& molecule1 = *it1;
				for (SingleCellIterator it2 = begin2; it2 != end2; ++it2) {
					Molecule& molecule2 = *it2;

					double dummy[3];
					double dd = molecule2.dist2(molecule1, dummy);

					if (dd < _cutoffRadiusSquare) {
	                    _rdf->observeRDF(molecule1, molecule2, dd);
					}
				}
			}
		} // isBoundaryCell
	}
}
