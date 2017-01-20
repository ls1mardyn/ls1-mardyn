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

	int numMolecules = getMoleculeCount();

	for (int i = 0; i < numMolecules; ++i) {
		if (moleculesAt(i).id() == molid) {
			found = true;
			deleteMoleculeByIndex(i);
			break;
		}
	}
	return found;
}
