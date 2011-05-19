/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "BlockTraverse.h"
#include "molecules/Molecule.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "Cell.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"

#include <vector>
#include <cmath>

using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################


BlockTraverse::BlockTraverse(
		ParticleContainer* moleculeContainer,
		vector<Cell>& cells, 
        vector<unsigned long>& innerCellIndices, 
        vector<unsigned long>& boundaryCellIndices,
        vector<unsigned long>& haloCellIndices,
		vector<vector<unsigned long> >& forwardNeighbourOffsets,
        vector<vector<unsigned long> >& backwardNeighbourOffsets
)
		: _moleculeContainer(moleculeContainer),
			_particlePairsHandler(moleculeContainer->getPairHandler()),
			_cells(cells),
			_innerCellIndices(innerCellIndices),
			_boundaryCellIndices(boundaryCellIndices),
            _haloCellIndices(haloCellIndices),
			_forwardNeighbourOffsets(&forwardNeighbourOffsets), _backwardNeighbourOffsets(&backwardNeighbourOffsets),
			_allocatedOffsets(false) {
}

BlockTraverse::BlockTraverse(
		ParticleContainer* moleculeContainer,
		vector<Cell>& cells, 
        vector<unsigned long>& innerCellIndices, 
        vector<unsigned long>& boundaryCellIndices, 
        vector<unsigned long>& haloCellIndices
)
		: _moleculeContainer(moleculeContainer),
			_particlePairsHandler(moleculeContainer->getPairHandler()),
			_cells(cells),
			_innerCellIndices(innerCellIndices),
			_boundaryCellIndices(boundaryCellIndices),
			_haloCellIndices(haloCellIndices),
			_forwardNeighbourOffsets(0), _backwardNeighbourOffsets(0),
			_allocatedOffsets(true)
{
	_forwardNeighbourOffsets = new vector<vector<unsigned long> >;
	_backwardNeighbourOffsets = new vector<vector<unsigned long> >;
}

BlockTraverse::~BlockTraverse() {
	if (_allocatedOffsets) {
		delete _forwardNeighbourOffsets;
		delete _backwardNeighbourOffsets;
	}
}

void BlockTraverse::assignOffsets(vector<unsigned long>& forwardNeighbourOffsets, vector<unsigned long>& backwardNeighbourOffsets) {
	_forwardNeighbourOffsets->assign(_cells.size(), forwardNeighbourOffsets);
	_backwardNeighbourOffsets->assign(_cells.size(), backwardNeighbourOffsets);
}

void BlockTraverse::traversePairs() {

	_particlePairsHandler = _moleculeContainer->getPairHandler();
	double _cutoffRadius = _moleculeContainer->getCutoff();
	double _LJCutoffRadius = _moleculeContainer->getLJCutoff();
	double _tersoffCutoffRadius = _moleculeContainer->getTersoffCutoff();
	vector<vector<unsigned long> >& forwardNeighbourOffsets = *_forwardNeighbourOffsets;
	vector<vector<unsigned long> >& backwardNeighbourOffsets = *_backwardNeighbourOffsets;

	_particlePairsHandler->init();

	// XXX comment
	double distanceVector[3];
	// loop over all cells
	vector<Molecule*>::iterator molIter1;
	vector<Molecule*>::iterator molIter2;

#ifndef NDEBUG
	// reset forces and momenta to zero
	global_log->debug() << "Resetting forces and momenta, disconnecting Tersoff pairs." << endl;
#endif
	{
		double zeroVec[3] = {0.0, 0.0, 0.0};

        // TODO: check if the reset is done twice as leaving this part has no difference on the result.
        Molecule *moleculePtr;
        for ( moleculePtr = _moleculeContainer->begin(); moleculePtr != _moleculeContainer->end(); moleculePtr = _moleculeContainer->next() ) {
				Molecule& molecule1 = *moleculePtr;
				molecule1.setF(zeroVec);
				molecule1.setM(zeroVec);
				molecule1.clearTersoffNeighbourList();
        }
	}

	vector<unsigned long>::iterator cellIndexIter;
	vector<unsigned long>::iterator neighbourOffsetsIter;

	// sqare of the cutoff radius
	double cutoffRadiusSquare = _cutoffRadius * _cutoffRadius;
	double LJCutoffRadiusSquare = _LJCutoffRadius * _LJCutoffRadius;
	double tersoffCutoffRadiusSquare = _tersoffCutoffRadius * _tersoffCutoffRadius;

#ifndef NDEBUG
	global_log->debug() << "Processing pairs and preprocessing Tersoff pairs." << endl;
#endif

	// loop over all inner cells and calculate forces to forward neighbours
	for (cellIndexIter = _innerCellIndices.begin(); cellIndexIter != _innerCellIndices.end(); cellIndexIter++) {
		unsigned long cellIndex = *cellIndexIter;
		Cell& currentCell = _cells[cellIndex];
        if( currentCell.getMoleculeCount() == 0 ) {
            continue;
        }

		// forces between molecules in the cell
		for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
			Molecule& molecule1 = **molIter1;
			unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching
			molIter2 = molIter1;
			molIter2++; // no self interaction

			for (; molIter2 != currentCell.getParticlePointers().end(); molIter2++) {
				Molecule& molecule2 = **molIter2;
				assert(&molecule1 != &molecule2);
				double dd = molecule2.dist2(molecule1, distanceVector);

				if (dd < cutoffRadiusSquare) {
					_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
					if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
					}
				}
			}
		}

		// loop over all neighbours
		for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
			Cell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];

            if( neighbourCell.getMoleculeCount() == 0 ) {
                continue;
            }
			// loop over all particles in the cell
			for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
				Molecule& molecule1 = **molIter1;
				unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

				for (molIter2 = neighbourCell.getParticlePointers().begin(); molIter2 != neighbourCell.getParticlePointers().end(); molIter2++) {
					Molecule& molecule2 = **molIter2;
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < cutoffRadiusSquare) {
						_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
							_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
						}
					}
				}

			}
		}
	}

	// loop over halo cells and detect Tersoff neighbours within the halo
	// this is relevant for the angle summation
    for (cellIndexIter = _haloCellIndices.begin(); cellIndexIter != _haloCellIndices.end(); cellIndexIter++) {
        unsigned long cellIndex = *cellIndexIter;
        Cell& currentCell = _cells[cellIndex];
        if( currentCell.getMoleculeCount() == 0 ) {
            continue;
        }

		for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
			Molecule& molecule1 = **molIter1;
			if (molecule1.numTersoff() == 0)
				continue;
			molIter2 = molIter1;
			molIter2++;
			for (; molIter2 != currentCell.getParticlePointers().end(); molIter2++) {
				Molecule& molecule2 = **molIter2;
				assert(&molecule1 != &molecule2);
				if (molecule2.numTersoff() > 0) {
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < tersoffCutoffRadiusSquare)
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, true);
				}
			}

			for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				int j = cellIndex + *neighbourOffsetsIter;
				if ((j < 0) || (j >= (int) (_cells.size())))
					continue;
				Cell& neighbourCell = _cells[j];
				if (!neighbourCell.isHaloCell())
					continue;
				for (molIter2 = neighbourCell.getParticlePointers().begin(); molIter2 != neighbourCell.getParticlePointers().end(); molIter2++) {
					Molecule& molecule2 = **molIter2;
					if (molecule2.numTersoff() == 0)
						continue;
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < tersoffCutoffRadiusSquare)
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, true);
				}
			}
		}
	}

	// loop over all boundary cells and calculate forces to forward and backward neighbours
	for (cellIndexIter = _boundaryCellIndices.begin(); cellIndexIter != _boundaryCellIndices.end(); cellIndexIter++) {
		unsigned long cellIndex = *cellIndexIter;
		Cell& currentCell = _cells[cellIndex];
        if( currentCell.getMoleculeCount() == 0 ) {
            continue;
        }

		// forces between molecules in the cell
		for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
			Molecule& molecule1 = **molIter1;
			unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching
			molIter2 = molIter1;
			molIter2++;
			for (; molIter2 != currentCell.getParticlePointers().end(); molIter2++) {
				Molecule& molecule2 = **molIter2;
				assert(&molecule1 != &molecule2);

				double dd = molecule2.dist2(molecule1, distanceVector);
				if (dd < cutoffRadiusSquare) {
					_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < LJCutoffRadiusSquare));
					if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
						_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, false);
					}
				}
			}
		}

		// loop over all forward neighbours
		for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
			Cell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];
            if( neighbourCell.getMoleculeCount() == 0 ) {
                continue;
            }

			// loop over all particles in the cell
			for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
				Molecule& molecule1 = **molIter1;
				unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

				for (molIter2 = neighbourCell.getParticlePointers().begin(); molIter2 != neighbourCell.getParticlePointers().end(); molIter2++) {
					Molecule& molecule2 = **molIter2;

					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < cutoffRadiusSquare) {
						PairType pairType = MOLECULE_MOLECULE;
						if (neighbourCell.isHaloCell() && ! molecule1.isLessThan(molecule2)) {
							/* Do not sum up values twice. */
							pairType = MOLECULE_HALOMOLECULE;
						}
						_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < LJCutoffRadiusSquare));
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
							_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
						}
					}
				}
			}
		}

		// loop over all backward neighbours. calculate only forces
		// to neighbour cells in the halo region, all others already have been calculated
		for (neighbourOffsetsIter = backwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter != backwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
			Cell& neighbourCell = _cells[cellIndex + *neighbourOffsetsIter];

			if (neighbourCell.isHaloCell() && neighbourCell.getMoleculeCount() > 0) {
				// loop over all particles in the cell
				for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
					Molecule& molecule1 = **molIter1;
					unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

					for (molIter2 = neighbourCell.getParticlePointers().begin(); molIter2 != neighbourCell.getParticlePointers().end(); molIter2++) {
						Molecule& molecule2 = **molIter2;

						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < cutoffRadiusSquare) {
							PairType pairType = molecule1.isLessThan(molecule2) ? MOLECULE_MOLECULE : MOLECULE_HALOMOLECULE;
							_particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < LJCutoffRadiusSquare));
							if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd < tersoffCutoffRadiusSquare)) {
								_particlePairsHandler->preprocessTersoffPair(molecule1, molecule2, (pairType == MOLECULE_HALOMOLECULE));
							}
						}
					}
				}
			}
		}
	}

#ifndef NDEBUG
	global_log->debug() << "processing Tersoff potential." << endl;
#endif
	double params[15];
	double delta_r = 0.;
	bool knowparams = false;

	for (cellIndexIter = _innerCellIndices.begin(); cellIndexIter != _boundaryCellIndices.end(); cellIndexIter++) {

		if (cellIndexIter == _innerCellIndices.end())
			cellIndexIter = _boundaryCellIndices.begin();

		Cell& currentCell = _cells[*cellIndexIter];
		for (molIter1 = currentCell.getParticlePointers().begin(); molIter1 != currentCell.getParticlePointers().end(); molIter1++) {
			Molecule& molecule1 = **molIter1;

			if (molecule1.numTersoff() == 0)
				continue;

			if (!knowparams) {
				delta_r = molecule1.tersoffParameters(params);
				knowparams = true;
			}
			_particlePairsHandler->processTersoffAtom(molecule1, params, delta_r);
		}
	}

	_particlePairsHandler->finish();
}
