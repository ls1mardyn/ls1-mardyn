//
// Created by alex on 31.01.24.
//

#include "TracerCellProcessor.h"
#include "particleContainer/ParticleCell.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "molecules/Molecule.h"

void TracerCellProcessor::processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) {
    double distanceVector[3];

    auto begin1 = cell1.iterator();
    auto begin2 = cell2.iterator();

    if(sumAll) { // sumAll - moleculesAt is gone, use SingleCellIterator now ?

        // loop over all particles in the cell
        for (auto it1 = begin1; it1.isValid(); ++it1) {
            Molecule& molecule1 = *it1;
            if(_comp_to_res[molecule1.componentid()] == CoarseGrain) continue;
            for (auto it2 = begin2; it2.isValid(); ++it2) {
                Molecule& molecule2 = *it2;
                if(_comp_to_res[molecule2.componentid()] == CoarseGrain) continue;
                if(molecule1.getID() == molecule2.getID()) continue;  // for grand canonical ensemble and traversal of pseudocells
                double dd = molecule2.dist2(molecule1, distanceVector);
                if (dd < _cutoffRadiusSquare) {
                    _particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < _LJCutoffRadiusSquare));
                }
            }

        }
    } else { // sumHalf
        if (cell1.isInnerCell()) {//no cell is halo
            // loop over all particles in the cell


            for (auto it1 = begin1; it1.isValid(); ++it1) {
                Molecule& molecule1 = *it1;
                if(_comp_to_res[molecule1.componentid()] == CoarseGrain) continue;
                for (auto it2 = begin2; it2.isValid(); ++it2) {
                    Molecule& molecule2 = *it2;
                    if(_comp_to_res[molecule2.componentid()] == CoarseGrain) continue;
                    if(molecule1.getID() == molecule2.getID()) continue;  // for grand canonical ensemble and traversal of pseudocells
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

            for (auto it1 = begin1; it1.isValid(); ++it1) {
                Molecule& molecule1 = *it1;
                if(_comp_to_res[molecule1.componentid()] == CoarseGrain) continue;
                for (auto it2 = begin2; it2.isValid(); ++it2) {
                    Molecule& molecule2 = *it2;
                    if(_comp_to_res[molecule2.componentid()] == CoarseGrain) continue;
                    double dd = molecule2.dist2(molecule1, distanceVector);
                    if (dd < _cutoffRadiusSquare) {
                        _particlePairsHandler->processPair(molecule1, molecule2, distanceVector, pairType, dd, (dd < _LJCutoffRadiusSquare));
                    }
                }
            }
        } // isBoundaryCell
    }
}

double TracerCellProcessor::processSingleMolecule(Molecule *m1, ParticleCell &cell2) {
    if(_comp_to_res[m1->componentid()] == CoarseGrain) return 0.0;

    double distanceVector[3];

    int neighbourParticleCount = cell2.getMoleculeCount();
    double u = 0.0;

    auto begin2 = cell2.iterator();

    for (auto it2 = begin2; it2.isValid(); ++it2) {
        Molecule& molecule2 = *it2;
        if(m1->getID() == molecule2.getID()) continue;
        if(_comp_to_res[molecule2.componentid()] == CoarseGrain) continue;

        double dd = molecule2.dist2(*m1, distanceVector);
        if (dd < _cutoffRadiusSquare)
        {
            PairType pairType = MOLECULE_MOLECULE_FLUID;
            u += _particlePairsHandler->processPair(*m1, molecule2, distanceVector, pairType, dd, (dd < _LJCutoffRadiusSquare));
        }
    }
    return u;
}

void TracerCellProcessor::processCell(ParticleCell &cell) {
    double distanceVector[3];

    if (cell.isInnerCell() || cell.isBoundaryCell()) {
        auto begin = cell.iterator();

        for (auto it1 = begin; it1.isValid(); ++it1) {
            Molecule& molecule1 = *it1;
            if(_comp_to_res[molecule1.componentid()] == CoarseGrain) continue;
            auto it2 = it1;
            ++it2;
            for (; it2.isValid(); ++it2) {
                Molecule& molecule2 = *it2;
                if(_comp_to_res[molecule2.componentid()] == CoarseGrain) continue;
                mardyn_assert(&molecule1 != &molecule2);
                double dd = molecule2.dist2(molecule1, distanceVector);

                if (dd < _cutoffRadiusSquare) {
                    _particlePairsHandler->processPair(molecule1, molecule2, distanceVector, MOLECULE_MOLECULE, dd, (dd < _LJCutoffRadiusSquare));
                }
            }
        }
    } // if (isInnerCell())
}

TracerCellProcessor::TracerCellProcessor(const double cutoffRadius, const double ljCutoffRadius,
                                         ParticlePairsHandler *particlePairsHandler, std::vector<Resolution>& ctr) : LegacyCellProcessor(cutoffRadius,
                                                                                                           ljCutoffRadius,
                                                                                                           particlePairsHandler),
                                                                                                         _comp_to_res(ctr), _particlePairsHandler(particlePairsHandler) {

}
