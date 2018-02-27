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
		auto begin = cell.iterator();

		for (auto it1 = begin; it1.hasNext(); it1.next()) {
			Molecule& molecule1 = *it1;

			auto it2 = it1;
			it2.next();
			for (; it2.hasNext(); it2.next()) {
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
	auto begin1 = cell1.iterator();
	auto begin2 = cell2.iterator();

	if(sumAll) { // sumAll - moleculesAt is gone, use auto now ?

		// loop over all particles in the cell
		for (auto it1 = begin1; it1.hasNext(); it1.next()) {
			Molecule& molecule1 = *it1; 
			for (auto it2 = begin2; it2.hasNext(); it2.next()) {
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


			for (auto it1 = begin1; it1.hasNext(); it1.next()) {
				Molecule& molecule1 = *it1;

				for (auto it2 = begin2; it2.hasNext(); it2.next()) {
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

			for (auto it1 = begin1; it1.hasNext(); it1.next()) {
				Molecule& molecule1 = *it1;
				for (auto it2 = begin2; it2.hasNext(); it2.next()) {
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
