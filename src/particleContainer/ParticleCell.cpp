#include "particleContainer/ParticleCell.h"
#include "molecules/Molecule.h"
#include "utils/UnorderedVector.h"

#include <cassert>
#include <vector>



using namespace std;

ParticleCell::ParticleCell() : _molecules(), _cellDataSoA(0) {
	for(int d = 0; d < 3; ++d) {
		_boxMin[d] = 0.0;
		_boxMax[d] = 0.0;
	}
}

ParticleCell::~ParticleCell() {
	assert(!_cellDataSoA);

	if(!isEmpty()) {
		deallocateAllParticles();
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

bool ParticleCell::addParticle(Molecule* particle_ptr, bool checkWhetherDuplicate) {
	bool wasInserted = false;

#ifndef NDEBUG
	bool isIn = particle_ptr->inBox(_boxMin, _boxMax);
//	assert(isIn);
#endif
	if (checkWhetherDuplicate == false) {
		_molecules.push_back(particle_ptr);
		wasInserted = true;
	} else {

		// perform a check whether this molecule exists (has been received) already

		vector<Molecule*>::const_iterator it;
		unsigned long pID = particle_ptr->id();

		bool found = false;

		for (it = _molecules.begin(); it != _molecules.end(); ++it) {
			if (pID == (*it)->id()) {
				found = true;
				break;
			}
		}

		if (not found) {
			_molecules.push_back(particle_ptr);
			wasInserted = true;
		}
	}
	return wasInserted;
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
			UnorderedVector::fastRemove(_molecules, molecule_iter);
			break;
		}
	}
	return found;
}

void ParticleCell::fastRemoveMolecule(std::vector<Molecule *>::iterator& it) {
	UnorderedVector::fastRemove(_molecules, it);
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
