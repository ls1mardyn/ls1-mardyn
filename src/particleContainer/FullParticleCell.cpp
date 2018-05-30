/*
 * FullParticleCell.cpp
 *
 *  Created on: 6 Feb 2017
 *      Author: tchipevn
 */

#include "particleContainer/ParticleCell.h"
#include "particleContainer/FullParticleCell.h"
#include "molecules/Molecule.h"
#include "utils/UnorderedVector.h"
#include "Simulation.h"

#include "utils/mardyn_assert.h"
#include <vector>

using namespace std;

FullParticleCell::FullParticleCell() :
		_molecules(), _cellDataSoA(0, 0, 0, 0, 0) { }

FullParticleCell::~FullParticleCell() {
//	if(!isEmpty()) {
//		deallocateAllParticles();
//	}
}

void FullParticleCell::deallocateAllParticles() {
	_molecules.clear();
}

bool FullParticleCell::addParticle(Molecule& particle, bool checkWhetherDuplicate) {
	bool wasInserted = false;

#ifndef NDEBUG
	bool isIn = particle.inBox(_boxMin, _boxMax);
	mardyn_assert(isIn);
#endif

	if (checkWhetherDuplicate == false) {
		_molecules.push_back(particle);
		wasInserted = true;
	} else {

		// perform a check whether this molecule exists (has been received) already

		bool found;
		size_t index;
		findMoleculeByID(found, index, particle.id());

		if (not found) {
			_molecules.push_back(particle);
			wasInserted = true;
		}
	}
	return wasInserted;
}

bool FullParticleCell::isEmpty() const {
	return _molecules.empty();
}

int FullParticleCell::getMoleculeCount() const {
	return _molecules.size();
}

bool FullParticleCell::deleteMoleculeByIndex(size_t index) {
//	mardyn_assert(index >= 0); - this is always true now
	mardyn_assert(index < _molecules.size());

	bool found = true;
	auto it = _molecules.begin() + index;
	UnorderedVector::fastRemove(_molecules, it);
	return found;
}

void FullParticleCell::preUpdateLeavingMolecules() {
	_leavingMolecules.clear();

	#ifndef NDEBUG
	const size_t size_total = _molecules.size(); // for debugging, see below
	#endif

	for (auto mol = _molecules.begin(); mol != _molecules.end(); ) {
		mol->setSoA(nullptr);

		const bool isStaying = mol->inBox(_boxMin, _boxMax);

		if (isStaying) {
			// don't do anything, just advance iterator
			++mol;
		}
		else {
			_leavingMolecules.push_back(*mol);
			UnorderedVector::fastRemove(_molecules, mol);
			// don't advance iterator, it now points to a new molecule
		}
	}

	mardyn_assert(_molecules.size() + _leavingMolecules.size() == size_total); // any molecules lost?
}

void FullParticleCell::updateLeavingMoleculesBase(ParticleCellBase& otherCell) {
	FullParticleCell& oCell = downcastReferenceFull(otherCell);
	updateLeavingMolecules(oCell);
}

void FullParticleCell::updateLeavingMolecules(FullParticleCell& otherCell){
	for (auto m = otherCell._leavingMolecules.begin(); m != otherCell._leavingMolecules.end(); ++m) { // loop over all indices
		if (m->inBox(_boxMin, _boxMax)) { // if molecule moves in this cell
			addParticle(*m);
		}
	}
}

void FullParticleCell::postUpdateLeavingMolecules(){
	_leavingMolecules.clear();
}

/*
std::vector<Molecule>& FullParticleCell::filterLeavingMolecules() {

	for (auto it = _molecules.begin(); it != _molecules.end();) {
		it->setSoA(nullptr);

		bool isStaying = it->inBox(_boxMin, _boxMax);

		if (isStaying) {
			++it; // next iteration
		} else {
			_leavingMolecules.push_back(*it);
			UnorderedVector::fastRemove(_molecules, it);
			// don't advance iterator, it now points to a new molecule
		}
	}

	return _leavingMolecules;
}
*/

void FullParticleCell::getRegion(double lowCorner[3], double highCorner[3], std::vector<Molecule*> &particlePtrs, bool removeFromContainer) {
	for (auto particleIter = _molecules.begin(); particleIter != _molecules.end();) {
		if (particleIter->inBox(lowCorner, highCorner)) {
			if (not removeFromContainer) {
				particlePtrs.push_back(&(*particleIter));
			} else {
				particlePtrs.push_back(new Molecule(*particleIter));
				UnorderedVector::fastRemove(_molecules, particleIter);
				// particleIter already points at next molecule, so continue without incrementing
				continue;
			}
		}
		++particleIter;
	}
}

void FullParticleCell::buildSoACaches() {

	// Determine the total number of centers.
	size_t numMolecules = _molecules.size();
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	size_t nDipoles = 0;
	size_t nQuadrupoles = 0;

	for (size_t m = 0;  m < numMolecules; ++m) {
		nLJCenters += _molecules[m].numLJcenters();
		nCharges += _molecules[m].numCharges();
		nDipoles += _molecules[m].numDipoles();
		nQuadrupoles += _molecules[m].numQuadrupoles();
	}

	// Construct the SoA.
	_cellDataSoA.resize(numMolecules,nLJCenters,nCharges,nDipoles,nQuadrupoles);

	size_t iLJCenters = 0;
	size_t iCharges = 0;
	size_t iDipoles = 0;
	size_t iQuadrupoles = 0;

	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < _molecules.size(); ++i) {
		Molecule & M = _molecules[i];
		const size_t mol_ljc_num = M.numLJcenters();
		const size_t mol_charges_num = M.numCharges();
		const size_t mol_dipoles_num = M.numDipoles();
		const size_t mol_quadrupoles_num = M.numQuadrupoles();

		_cellDataSoA._mol_ljc_num[i] = mol_ljc_num;
		_cellDataSoA._mol_charges_num[i] = mol_charges_num;
		_cellDataSoA._mol_dipoles_num[i] = mol_dipoles_num;
		_cellDataSoA._mol_quadrupoles_num[i] = mol_quadrupoles_num;

		_cellDataSoA._mol_pos.x(i) = M.r(0);
		_cellDataSoA._mol_pos.y(i) = M.r(1);
		_cellDataSoA._mol_pos.z(i) = M.r(2);

		M.setupSoACache(&_cellDataSoA, iLJCenters, iCharges, iDipoles, iQuadrupoles);

		iLJCenters += mol_ljc_num;
		iCharges += mol_charges_num;
		iDipoles += mol_dipoles_num;
		iQuadrupoles += mol_quadrupoles_num;

		M.clearFM();
	}
}

void FullParticleCell::reserveMoleculeStorage(size_t numMols) {
	_molecules.reserve(numMols);
}
