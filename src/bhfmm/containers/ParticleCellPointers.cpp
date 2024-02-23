/*
 * ParticleCellPointers.cpp
 *
 *  Created on: 20 Sep 2016
 *      Author: tchipevn
 */

#include "ParticleCellPointers.h"
#include "molecules/Molecule.h"
#include "utils/UnorderedVector.h"

#include "utils/mardyn_assert.h"
#include <vector>

namespace bhfmm {

ParticleCellPointers::ParticleCellPointers() :
		haloCell(false), boundaryCell(false), innerCell(false), innerMostCell(false),
		_molecules(), _cellDataSoA(0, 0, 0, 0, 0) {
	// TODO: disable when in MoleculeRMM mode. It's only 1CLJ anyway.
	for (int d = 0; d < 3; ++d) {
		_boxMin[d] = 0.0;
		_boxMax[d] = 0.0;
	}
}

ParticleCellPointers::~ParticleCellPointers() {
//	if(!isEmpty()) {
//		deallocateAllParticles();
//	}
}

void ParticleCellPointers::removeAllParticles() {
	_molecules.clear();
}

bool ParticleCellPointers::addParticle(Molecule* particle_ptr) {
#ifndef NDEBUG
	bool isIn = particle_ptr->inBox(_boxMin, _boxMax);
	mardyn_assert(isIn);
#endif

	_molecules.push_back(particle_ptr);
	bool wasInserted = true;

	return wasInserted;
}

bool ParticleCellPointers::isEmpty() const {
	return _molecules.empty();
}

int ParticleCellPointers::getMoleculeCount() const {
	return _molecules.size();
}

} /* namespace bhfmm */
