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

bool ParticleContainer::isInBoundingBox(double r[3]) const {
	if (r[0] >= _boundingBoxMin[0] && r[1] >= _boundingBoxMin[1] && r[2] >= _boundingBoxMin[2]
			&& r[0] < _boundingBoxMax[0] && r[1] < _boundingBoxMax[1] && r[2] < _boundingBoxMax[2]) {
		return true;
	} else {
		return false;
	}
}

int ParticleContainer::getHaloWidthNumCells() {
	return 0;
}
