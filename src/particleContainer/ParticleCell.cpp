#include "particleContainer/ParticleCell.h"

#include <cassert>
#include <vector>

#include "molecules/Molecule.h"

using namespace std;

ParticleCell::ParticleCell() : _molecules(), _cellDataSoA(0) {
	for(int d = 0; d < 3; ++d) {
		_boxMin[d] = 0.0;
		_boxMax[d] = 0.0;
	}
}

ParticleCell::~ParticleCell() {
	assert(!_cellDataSoA);

	// ugly, but necessary, because the AdaptiveSubCells still work with a global list of references
	if(!isEmpty()) {
		removeAllParticles();
	}
}

void ParticleCell::removeAllParticles() {
	_molecules.clear();
}

void ParticleCell::deallocateAllParticles() {
	std::vector<Molecule * >::iterator it;
	for(it = _molecules.begin(); it != _molecules.end(); ++it) {
		delete *it;
	}
	removeAllParticles();
}

void ParticleCell::addParticle(Molecule* particle_ptr) {
#ifndef NDEBUG
	bool isIn = particle_ptr->inBox(_boxMin, _boxMax);
//	assert(isIn);
#endif
	_molecules.push_back(particle_ptr);
}

vector<Molecule*>& ParticleCell::getParticlePointers() {
	return _molecules;
}

const vector<Molecule*>& ParticleCell::getParticlePointers() const {
	return _molecules;
}

bool ParticleCell::isEmpty() const {
	return _molecules.empty();
}

int ParticleCell::getMoleculeCount() const {
	return _molecules.size();
}

bool ParticleCell::deleteMolecule(unsigned long molecule_id) {
	bool found = false;
	vector<Molecule*>::iterator molecule_iter;

	for (molecule_iter = _molecules.begin(); molecule_iter != _molecules.end(); molecule_iter++) {
		Molecule *molecule = *molecule_iter;
		if (molecule->id() == molecule_id) {
			found = true;
			_molecules.erase(molecule_iter);
			break;
		}
	}
	return found;
}

std::vector<Molecule*>& ParticleCell::filterLeavingMolecules() {
	std::vector<Molecule*>::iterator it;

	for (it = _molecules.begin(); it != _molecules.end();) {
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
