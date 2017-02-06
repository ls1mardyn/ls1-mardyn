#ifndef PARTICLE_CELL_H_
#define PARTICLE_CELL_H_

/**
 * the old class ParticleCell is now called FullParticleCell and it implements ParticleCellBase.
 */

#include "FullParticleCell.h"
#include "ParticleCellWR.h"

#ifndef MARDYN_WR
	typedef FullParticleCell ParticleCell;
#else
	typedef ParticleCell_WR ParticleCell;
#endif

inline FullParticleCell* downcastPointerFull(ParticleCellBase* c) {
	assert(static_cast<FullParticleCell*>(c) == dynamic_cast<FullParticleCell*>(c));
	return static_cast<FullParticleCell*>(c);
}

inline FullParticleCell& downcastReferenceFull(ParticleCellBase& c) {
	assert(&static_cast<FullParticleCell&>(c) == &dynamic_cast<FullParticleCell&>(c));
	return static_cast<FullParticleCell&>(c);
}

#endif /* PARTICLE_CELL_H_ */
