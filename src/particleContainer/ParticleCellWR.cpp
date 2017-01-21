/*
 * ParticleCellWR.cpp
 *
 *  Created on: 20 Jan 2017
 *      Author: tchipevn
 */

#include "ParticleCellWR.h"

ParticleCell_WR::ParticleCell_WR() : _cellDataSoA_WR(0) {
	// TODO Auto-generated constructor stub

}

ParticleCell_WR::~ParticleCell_WR() {
	// TODO Auto-generated destructor stub
}

void ParticleCell_WR::deallocateAllParticles() {
	_cellDataSoA_WR.resize(0);
}

bool ParticleCell_WR::addParticle(Molecule& particle, bool checkWhetherDuplicate) {

	bool wasInserted;
	bool found = false;

	if (checkWhetherDuplicate) {
		size_t index;
		findMoleculeByID(found, index, particle.id());
	}

	if (not found) {
		_cellDataSoA_WR.appendMolecule(particle);
		wasInserted = true;
	} else {
		wasInserted = false;
	}

	return wasInserted;
}

Molecule& ParticleCell_WR::moleculesAt(size_t i) {
	_cellDataSoA_WR.readMolecule(i, _dummy);
	return _dummy;
}

bool ParticleCell_WR::isEmpty() const {
	return getMoleculeCount() == 0;
}

bool ParticleCell_WR::deleteMoleculeByIndex(size_t index) {
	assert(false);
	return false;
}

int ParticleCell_WR::getMoleculeCount() const {
	return _cellDataSoA_WR._mol_num;
}

void ParticleCell_WR::preUpdateLeavingMolecules() {
}

void ParticleCell_WR::updateLeavingMoleculesBase(ParticleCellBase& otherCell) {
}

void ParticleCell_WR::postUpdateLeavingMolecules() {
}

void ParticleCell_WR::getRegion(double lowCorner[3], double highCorner[3],
		std::vector<Molecule*>& particlePtrs, bool removeFromContainer) {
}

void ParticleCell_WR::buildSoACaches() {
}

void ParticleCell_WR::reserveMoleculeStorage(size_t numMols) {
	_cellDataSoA_WR.resize(numMols);
}
