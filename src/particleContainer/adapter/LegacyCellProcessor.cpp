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

	for (int j = 0; j < neighbourParticleCount; j++) {
		Molecule& molecule2 = cell2.moleculesAt(j);
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

int LegacyCellProcessor::countNeighbours(Molecule* m1, ParticleCell& cell2, double RR)
{
        int tn = 0;
        double distanceVector[3];

        int neighbourParticleCount = cell2.getMoleculeCount();

        for (int j = 0; j < neighbourParticleCount; j++) {
                Molecule& molecule2 = cell2.moleculesAt(j);
                if(m1->id() == molecule2.id()) continue;
                double dd = molecule2.dist2(*m1, distanceVector);
                if (dd < RR) tn++;
        }
        return tn;
}

void LegacyCellProcessor::processCellPair(ParticleCell& cell1, ParticleCell& cell2) {
	double distanceVector[3];

	int currentParticleCount = cell1.getMoleculeCount();
	int neighbourParticleCount = cell2.getMoleculeCount();

	if (cell1.isInnerCell()) {//no cell is halo
		// loop over all particles in the cell
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = cell1.moleculesAt(i);

			for (int j = 0; j < neighbourParticleCount; j++) {
				Molecule& molecule2 = cell2.moleculesAt(j);
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
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = cell1.moleculesAt(i);
			for (int j = 0; j < neighbourParticleCount; j++) {
				Molecule& molecule2 = cell2.moleculesAt(j);

				double dd = molecule2.dist2(molecule1, distanceVector);
				if (dd < _cutoffRadiusSquare) {
					PairType pairType = MOLECULE_MOLECULE;
					//if (cell2.isHaloCell() && ! molecule1.isLessThan(molecule2)) {//boundary <-> halo: using macroscopic boundary condition
					if (cell2.isHaloCell() && ! (cell1.getCellIndex()<cell2.getCellIndex())) {//boundary <-> halo: not using macroscopic boundary condition, instead cell indices are compared.
						/* Do not sum up values twice. */
						pairType = MOLECULE_HALOMOLECULE;
					}
					_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < _LJCutoffRadiusSquare));
				}
			}
		}
	} // isBoundaryCell
}

void LegacyCellProcessor::processCell(ParticleCell& cell) {
	double distanceVector[3];
	int currentParticleCount = cell.getMoleculeCount();

	if (cell.isInnerCell() || cell.isBoundaryCell()) {
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = cell.moleculesAt(i);

			for (int j = i+1; j < currentParticleCount; j++) {
				Molecule& molecule2 = cell.moleculesAt(j);
				assert(&molecule1 != &molecule2);
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
