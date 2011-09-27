#include <iostream>

#include "particleContainer/ParticleContainer.h"

#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "molecules/Molecule.h"


using namespace std;

ParticleContainer::ParticleContainer(ParticlePairsHandler* partPairsHandler, double bBoxMin[3], double bBoxMax[3])
		: _particlePairsHandler(partPairsHandler) {
	for (int i = 0; i < 3; i++) {
		_boundingBoxMin[i] = bBoxMin[i];
		_boundingBoxMax[i] = bBoxMax[i];
	}
}

ParticleContainer::~ParticleContainer() {
}

void ParticleContainer::rebuild(double bBoxMin[3], double bBoxMax[3]) {
	cout << "REBUILD OF PARTICLE CONTAINER" << endl;
	for (int i = 0; i < 3; i++) {
		_boundingBoxMin[i] = bBoxMin[i];
		_boundingBoxMax[i] = bBoxMax[i];
	}
}

double ParticleContainer::getBoundingBoxMin(int dimension) const {
	return this->_boundingBoxMin[dimension];
}

double ParticleContainer::getBoundingBoxMax(int dimension) const {
	return this->_boundingBoxMax[dimension];
}

double ParticleContainer::get_halo_L(int index) const {
	cerr << "ERROR: ParticleContainer::get_halo_L(...) has to be implemented in derived class" << endl;
	return 0;
}

void ParticleContainer::setPairHandler(ParticlePairsHandler* partPairHandler) {
	_particlePairsHandler = partPairHandler;
}

ParticlePairsHandler* ParticleContainer::getPairHandler() {
	return _particlePairsHandler;
}

void ParticleContainer::updateMoleculeCaches() {
	Molecule *tM;
	for (tM = this->begin(); tM != this->end(); tM = this->next() ) {
		tM->upd_cache();
	}
}
