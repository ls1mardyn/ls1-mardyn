/*
 * LegacyCellProcessor.cpp
 *
 * @Date: 18.03.2012
 * @Author: eckhardw
 */

#include "LegacyCellProcessor.h"
#include "particleContainer/ParticleCell.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include <vector>

using namespace std;
using namespace Log;

LegacyCellProcessor::LegacyCellProcessor(const double cutoffRadius, const double LJCutoffRadius,
		ParticlePairsHandler* particlePairsHandler)
: CellProcessor(cutoffRadius, LJCutoffRadius),
  _particlePairsHandler(particlePairsHandler){
}


LegacyCellProcessor::~LegacyCellProcessor() {
}


void LegacyCellProcessor::initTraversal()
{
	_particlePairsHandler->init();
}

double LegacyCellProcessor::processSingleMolecule(Molecule* m1, ParticleCell& cell2)
{
	double distanceVector[3];

	int neighbourParticleCount = cell2.getMoleculeCount();
	double u = 0.0;

	SingleCellIterator begin2 = cell2.iteratorBegin();
	SingleCellIterator end2 = cell2.iteratorEnd();

	for (SingleCellIterator it2 = begin2; it2 != end2; ++it2) {
		Molecule& molecule2 = *it2;
		if(m1->id() == molecule2.id()) continue;
		double dd = molecule2.dist2(*m1, distanceVector);
		if (dd < _cutoffRadiusSquare)
		{
			PairType pairType = MOLECULE_MOLECULE_FLUID;
			u += _particlePairsHandler->processPair(*m1, molecule2, distanceVector, pairType, dd, (dd < _LJCutoffRadiusSquare));
		}
	}
	return u;
}

void LegacyCellProcessor::processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll /* = false */) {
	double distanceVector[3];
	
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
				if(molecule1.id() == molecule2.id()) continue;  // for grand canonical ensemble and traversal of pseudocells
				double dd = molecule2.dist2(molecule1, distanceVector);
				if (dd < _cutoffRadiusSquare) {
					_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < _LJCutoffRadiusSquare));
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

					if(molecule1.id() == molecule2.id()) continue;  // for grand canonical ensemble and traversal of pseudocells
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < _cutoffRadiusSquare) {
						_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < _LJCutoffRadiusSquare));
					}
				}

			}
		} // inner cell

		if (cell1.isBoundaryCell()) {//first cell is boundary
			// loop over all particles in the cell
			PairType pairType = MOLECULE_MOLECULE;
			//if (cell2.isHaloCell() && ! molecule1.isLessThan(molecule2)) {//boundary <-> halo: using macroscopic boundary condition
			if (cell2.isHaloCell() && ! (cell1.getCellIndex()<cell2.getCellIndex())) {//boundary <-> halo: not using macroscopic boundary condition, instead cell indices are compared.
				/* Do not sum up values twice. */
				pairType = MOLECULE_HALOMOLECULE;
			}

			for (SingleCellIterator it1 = begin1; it1 != end1; ++it1) {
				Molecule& molecule1 = *it1;
				for (SingleCellIterator it2 = begin2; it2 != end2; ++it2) {
					Molecule& molecule2 = *it2;
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < _cutoffRadiusSquare) {
						_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < _LJCutoffRadiusSquare));
					}
				}
			}
		} // isBoundaryCell
	}
}



#if defined(ENABLE_REDUCED_MEMORY_MODE) && ENABLE_VECTORIZED_CODE==0
#warning Reduced memory mode (RMM) should now work using LegacyCellProcessor with the introduction of SingleCellIterator-s. Please double check at first usage!
#endif
////TODO: fix LegacyCellProcessor for RMM mode - especially the calls to moleculesAt(i) using the _dummy are dangerous!
// should work now, with SingleCellIterator

void LegacyCellProcessor::processCell(ParticleCell& cell) {
	double distanceVector[3];

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
				double dd = molecule2.dist2(molecule1, distanceVector);

				if (dd < _cutoffRadiusSquare) {
					_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < _LJCutoffRadiusSquare));
				}
			}
		}
	} // if (isInnerCell())
}

void LegacyCellProcessor::endTraversal() {
	_particlePairsHandler->finish();
}
