/*
 * ParticleCellBase.cpp
 *
 *  Created on: 20 Jan 2017
 *      Author: tchipevn
 */

#include "ParticleCellBase.h"

ParticleCellBase::ParticleCellBase() {
	// TODO Auto-generated constructor stub

}

ParticleCellBase::~ParticleCellBase() {
	// TODO Auto-generated destructor stub
}

bool ParticleCellBase::deleteMoleculeByID(unsigned long molid) {
	bool found = false;

	size_t index;
	findMoleculeByID(found, index, molid);
	if (found) {
		deleteMoleculeByIndex(index);
	}

	return found;
}

void ParticleCellBase::findMoleculeByID(bool& wasFound, size_t& index, unsigned long molid) const {
	wasFound = false;
	int numMolecules = getMoleculeCount();

	for (int i = 0; i < numMolecules; ++i) {
		if (moleculesAtConst(i).id() == molid) {
			index = i;
			wasFound = true;
			break;
		}
	}
}
