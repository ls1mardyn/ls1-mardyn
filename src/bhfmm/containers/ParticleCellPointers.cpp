/*
 * ParticleCellPointers.cpp
 *
 *  Created on: 20 Sep 2016
 *      Author: tchipevn
 */

#include "ParticleCellPointers.h"
#include "molecules/Molecule.h"
#include "utils/UnorderedVector.h"

#include <cassert>
#include <vector>

using namespace std;

namespace bhfmm {

ParticleCellPointers::ParticleCellPointers() :
		_molecules(), _cellDataSoA(0, 0, 0, 0, 0) { }

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
	assert(isIn);
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
