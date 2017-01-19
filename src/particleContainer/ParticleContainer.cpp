#include "particleContainer/ParticleContainer.h"

#include "molecules/Molecule.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

ParticleContainer::ParticleContainer(double bBoxMin[3], double bBoxMax[3]) {
	for (int i = 0; i < 3; i++) {
		_boundingBoxMin[i] = bBoxMin[i];
		_boundingBoxMax[i] = bBoxMax[i];
	}
}

ParticleContainer::~ParticleContainer() {
}

void ParticleContainer::rebuild(double bBoxMin[3], double bBoxMax[3]) {
	global_log->info() << "REBUILD OF PARTICLE CONTAINER" << endl;
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
double ParticleContainer::getHaloWidthNumCells() {
	return 0;
}
void ParticleContainer::updateMoleculeCaches() {
	Molecule *tM;
	for (tM = this->begin(); tM != this->end(); tM = this->next() ) {
		tM->upd_cache();
		tM->clearFM();
	}
}
