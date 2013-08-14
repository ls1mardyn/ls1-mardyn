#include "particleContainer/ParticleCell.h"

#include <cassert>
#include <vector>

#include "molecules/Molecule.h"

using namespace std;

ParticleCell::ParticleCell() :
		_cellDataSoA(0) {
}

ParticleCell::~ParticleCell() {
	assert(!_cellDataSoA);
}

void ParticleCell::removeAllParticles() {
	molecules.clear();
}

void ParticleCell::addParticle(Molecule* particle_ptr) {
	molecules.push_back(particle_ptr);
}

vector<Molecule*>& ParticleCell::getParticlePointers() {
	return molecules;
}

int ParticleCell::getMoleculeCount() const {
	return molecules.size();
}

bool ParticleCell::deleteMolecule(unsigned long molecule_id) {
	bool found = false;
	vector<Molecule*>::iterator molecule_iter;

	for (molecule_iter = molecules.begin(); molecule_iter != molecules.end(); molecule_iter++) {
		Molecule *molecule = *molecule_iter;
		if (molecule->id() == molecule_id) {
			found = true;
			molecules.erase(molecule_iter);
			break;
		}
	}
	return found;
}
