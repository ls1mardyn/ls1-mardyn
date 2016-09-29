#include "particleContainer/ParticleCell.h"
#include "molecules/Molecule.h"
#include "utils/UnorderedVector.h"

#include <cassert>
#include <vector>



using namespace std;

ParticleCell::ParticleCell() :
		_molecules(), _cellDataSoA(0, 0, 0, 0, 0), _cellIndex(0) {
	for (int d = 0; d < 3; ++d) {
		_boxMin[d] = 0.0;
		_boxMax[d] = 0.0;
	}
}

ParticleCell::~ParticleCell() {
//	if(!isEmpty()) {
//		deallocateAllParticles();
//	}
}

void ParticleCell::removeAllParticles() {
	_molecules.clear();
}

void ParticleCell::deallocateAllParticles() {
	for(auto it = _molecules.begin(); it != _molecules.end(); ++it) {
		delete *it;
	}
	removeAllParticles();
}

bool ParticleCell::addParticle(Molecule& particle, bool checkWhetherDuplicate) {
	bool wasInserted = false;

#ifndef NDEBUG
	bool isIn = particle.inBox(_boxMin, _boxMax);
	assert(isIn);
#endif

	if (checkWhetherDuplicate == false) {
		_molecules.push_back(new Molecule(particle));
		wasInserted = true;
	} else {

		// perform a check whether this molecule exists (has been received) already

		vector<Molecule*>::const_iterator it;
		unsigned long pID = particle.id();

		bool found = false;

		for (it = _molecules.begin(); it != _molecules.end(); ++it) {
			if (pID == (*it)->id()) {
				found = true;
				break;
			}
		}

		if (not found) {
			_molecules.push_back(new Molecule(particle));
			wasInserted = true;
		}
	}
	return wasInserted;
}

bool ParticleCell::isEmpty() const {
	return _molecules.empty();
}

int ParticleCell::getMoleculeCount() const {
	return _molecules.size();
}

bool ParticleCell::deleteMoleculeByID(unsigned long molecule_id) {
	bool found = false;
	vector<Molecule*>::iterator molecule_iter;

	for (molecule_iter = _molecules.begin(); molecule_iter != _molecules.end(); molecule_iter++) {
		Molecule *molecule = *molecule_iter;
		if (molecule->id() == molecule_id) {
			found = true;
			UnorderedVector::fastRemove(_molecules, molecule_iter);
			break;
		}
	}
	return found;
}

bool ParticleCell::deleteMoleculeByIndex(std::vector<Molecule * >::size_type index) {
//	assert(index >= 0); - this is always true now
	assert(index < _molecules.size());

	bool found = true;
	std::vector<Molecule *>::iterator it = _molecules.begin() + index;
	delete *it;
	UnorderedVector::fastRemove(_molecules, it);
	return found;
}

std::vector<Molecule*>& ParticleCell::filterLeavingMolecules() {
	std::vector<Molecule*>::iterator it;

	for (it = _molecules.begin(); it != _molecules.end();) {
		(*it)->setSoA(nullptr);
		const Molecule * const mol = *it;

		bool isStaying = mol->inBox(_boxMin, _boxMax);

		if (isStaying) {
			++it; // next iteration
		} else {
			_leavingMolecules.push_back(*it);

			// erase from _molecules:
			*it = _molecules.back(); // iterator points to a new molecule
			_molecules.back() = NULL;
			_molecules.pop_back();
			// don't advance iterator, it now points to a new molecule
		}
	}

	return _leavingMolecules;
}

void ParticleCell::getRegion(double lowCorner[3], double highCorner[3], std::vector<Molecule*> &particlePtrs, bool removeFromContainer) {
	std::vector<Molecule *>::iterator particleIter;

	for (particleIter = _molecules.begin(); particleIter != _molecules.end();) {
		if ((*particleIter)->inBox(lowCorner, highCorner)) {
			particlePtrs.push_back(*particleIter);
			if (removeFromContainer) {
				UnorderedVector::fastRemove(_molecules, particleIter);
				// particleIter already points at next molecule, so continue without incrementing
				continue;
			}
		}
		++particleIter;
	}
}

void ParticleCell::buildSoACaches() {

	// Determine the total number of centers.
	size_t numMolecules = _molecules.size();
	size_t nLJCenters = 0;
	size_t nCharges = 0;
	size_t nDipoles = 0;
	size_t nQuadrupoles = 0;

	for (size_t m = 0;  m < numMolecules; ++m) {
		nLJCenters += _molecules[m]->numLJcenters();
		nCharges += _molecules[m]->numCharges();
		nDipoles += _molecules[m]->numDipoles();
		nQuadrupoles += _molecules[m]->numQuadrupoles();
	}

	// Construct the SoA.
	_cellDataSoA.resize(numMolecules,nLJCenters,nCharges,nDipoles,nQuadrupoles);

	size_t iLJCenters = 0;
	size_t iCharges = 0;
	size_t iDipoles = 0;
	size_t iQuadrupoles = 0;

	// For each molecule iterate over all its centers.
	for (size_t i = 0; i < _molecules.size(); ++i) {
		Molecule & M = *_molecules[i];
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
