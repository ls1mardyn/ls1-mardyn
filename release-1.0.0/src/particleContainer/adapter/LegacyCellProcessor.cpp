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
		const double tersoffCutoffRadius, ParticlePairsHandler* particlePairsHandler)
: _cutoffRadiusSquare(cutoffRadius * cutoffRadius), _LJCutoffRadiusSquare(LJCutoffRadius * LJCutoffRadius),
  _tersoffCutoffRadiusSquare(tersoffCutoffRadius*tersoffCutoffRadius), _particlePairsHandler(particlePairsHandler){
	  /** @todo Check for multiple tersoff potentials with different parameters as the LegacyCellProcessor::postprocessCell() does only use one parameter set. */
	  global_log->warning() << "Note: The LegacyCellProcessor does not support multiple Tersoff sites with different parameters." << endl;
}


LegacyCellProcessor::~LegacyCellProcessor() {
}


void LegacyCellProcessor::initTraversal(const size_t numCells)
{
	_particlePairsHandler->init();
}


void LegacyCellProcessor::preprocessCell(ParticleCell& cell) {
	assert(!cell.isInActiveWindow());
	cell.setInActiveWindow();

	double zeroVec[3] = {0.0, 0.0, 0.0};

	// TODO: check if the reset is done twice as leaving this part has no difference on the result.
	vector<Molecule*>& particlePointers = cell.getParticlePointers();
	size_t size = cell.getParticlePointers().size();
	for (size_t i = 0; i < size; i++) {
		Molecule& molecule1 = *particlePointers[i];
		molecule1.setF(zeroVec);
		molecule1.setM(zeroVec);
		molecule1.clearTersoffNeighbourList();
	}
}

double LegacyCellProcessor::processSingleMolecule(Molecule* m1, ParticleCell& cell2)
{
	assert(cell2.isInActiveWindow());
	double distanceVector[3];

	std::vector<Molecule*>& neighbourCellParticles = cell2.getParticlePointers();
	int neighbourParticleCount = neighbourCellParticles.size();
	double u = 0.0;

	for (int j = 0; j < neighbourParticleCount; j++) {
		Molecule& molecule2 = *neighbourCellParticles[j];
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

void LegacyCellProcessor::processCellPair(ParticleCell& cell1, ParticleCell& cell2) {
	assert(cell1.isInActiveWindow());
	assert(cell2.isInActiveWindow());
	double distanceVector[3];

	std::vector<Molecule*>& currentCellParticles = cell1.getParticlePointers();
	int currentParticleCount = currentCellParticles.size();
	std::vector<Molecule*>& neighbourCellParticles = cell2.getParticlePointers();
	int neighbourParticleCount = neighbourCellParticles.size();

	if (cell1.isInnerCell()) {
		// loop over all particles in the cell
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = *currentCellParticles[i];
			unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

			for (int j = 0; j < neighbourParticleCount; j++) {
				Molecule& molecule2 = *neighbourCellParticles[j];
				double dd = molecule2.dist2(molecule1, distanceVector);
				if (dd < _cutoffRadiusSquare) {
					_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < _LJCutoffRadiusSquare));
					if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < _tersoffCutoffRadiusSquare)) {
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
					}
				}
			}

		}
	} // inner cell


	if (cell1.isHaloCell()) {
		assert(cell2.isHaloCell());
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = *currentCellParticles[i];
			if (molecule1.numTersoff() == 0) {
				continue;
			}
			for (int j = 0; j < neighbourParticleCount; j++) {
				Molecule& molecule2 = *neighbourCellParticles[j];
				assert(&molecule1 != &molecule2);
				if (molecule2.numTersoff() > 0) {
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < _tersoffCutoffRadiusSquare) {
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, true);
					}
				}
			}
		}
	} // isHaloCell


	if (cell1.isBoundaryCell()) {
		// loop over all particles in the cell
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = *currentCellParticles[i];
			unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

			for (int j = 0; j < neighbourParticleCount; j++) {
				Molecule& molecule2 = *neighbourCellParticles[j];

				double dd = molecule2.dist2(molecule1, distanceVector);
				if (dd < _cutoffRadiusSquare) {
					PairType pairType = MOLECULE_MOLECULE;
					if (cell2.isHaloCell() && ! molecule1.isLessThan(molecule2)) {
						/* Do not sum up values twice. */
						pairType = MOLECULE_HALOMOLECULE;
					}
					_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < _LJCutoffRadiusSquare));
					if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < _tersoffCutoffRadiusSquare)) {
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
					}
				}
			}
		}
	} // if ( isBoundaryCell() )
}

void LegacyCellProcessor::processCell(ParticleCell& cell) {
	assert(cell.isInActiveWindow());
	double distanceVector[3];
	std::vector<Molecule*>& currentCellParticles = cell.getParticlePointers();
	int currentParticleCount = currentCellParticles.size();

	if (cell.isInnerCell() || cell.isBoundaryCell()) {
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = *currentCellParticles[i];
			unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

			for (int j = i+1; j < currentParticleCount; j++) {
				Molecule& molecule2 = *currentCellParticles[j];
				assert(&molecule1 != &molecule2);
				double dd = molecule2.dist2(molecule1, distanceVector);

				if (dd < _cutoffRadiusSquare) {
					_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < _LJCutoffRadiusSquare));
					if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < _tersoffCutoffRadiusSquare)) {
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
					}
				}
			}
		}
	} // if (isInnerCell())

	if (cell.isHaloCell()) {
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = *currentCellParticles[i];
			if (molecule1.numTersoff() == 0)
				continue;

			for (int j = i+1; j < currentParticleCount; j++) {
				Molecule& molecule2 = *currentCellParticles[j];
				assert(&molecule1 != &molecule2);
				if (molecule2.numTersoff() > 0) {
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < _tersoffCutoffRadiusSquare)
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, true);
				}
			}
		}
	} // isHaloCell
}

void LegacyCellProcessor::postprocessCell(ParticleCell& cell) {
	assert(cell.isInActiveWindow());
	cell.clearInActiveWindow();

	std::vector<Molecule*>& currentCellParticles = cell.getParticlePointers();
	int currentParticleCount = currentCellParticles.size();

	double params[15];
	double delta_r = 0.;
	bool knowparams = false;

	if (cell.isInnerCell() || cell.isBoundaryCell()) {
		for (int i = 0; i < currentParticleCount; i++) {
			Molecule& molecule1 = *currentCellParticles[i];
			/** @todo Fix potential bug in case of different Tersoff sites with different parameters. */
			if (molecule1.numTersoff() > 0) {
				if (!knowparams) {
					delta_r = molecule1.tersoffParameters(params);
					knowparams = true;
				}
				_particlePairsHandler->processTersoffAtom(molecule1, params, delta_r);
			}
			molecule1.calcFM();
		}
	}
}

void LegacyCellProcessor::endTraversal() {
	_particlePairsHandler->finish();
}
