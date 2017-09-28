#ifndef PARTICLE_CELL_H_
#define PARTICLE_CELL_H_

/**
 * the old class ParticleCell is now called FullParticleCell and it implements ParticleCellBase.
 */

#include "FullParticleCell.h"
#include "ParticleCellWR.h"

#ifndef ENABLE_REDUCED_MEMORY_MODE
	typedef FullParticleCell ParticleCell;
#else
	typedef ParticleCell_WR ParticleCell;
#endif

inline FullParticleCell* downcastCellPointerFull(ParticleCellBase* c) {
	mardyn_assert(static_cast<FullParticleCell*>(c) == dynamic_cast<FullParticleCell*>(c));
	return static_cast<FullParticleCell*>(c);
}

inline FullParticleCell& downcastCellReferenceFull(ParticleCellBase& c) {
	mardyn_assert(&static_cast<FullParticleCell&>(c) == &dynamic_cast<FullParticleCell&>(c));
	return static_cast<FullParticleCell&>(c);
}

inline ParticleCell_WR* downcastCellPointerWR(ParticleCellBase* c) {
	mardyn_assert(static_cast<ParticleCell_WR*>(c) == dynamic_cast<ParticleCell_WR*>(c));
	return static_cast<ParticleCell_WR*>(c);
}

inline ParticleCell_WR& downcastCellReferenceWR(ParticleCellBase& c) {
	mardyn_assert(&static_cast<ParticleCell_WR&>(c) == &dynamic_cast<ParticleCell_WR&>(c));
	return static_cast<ParticleCell_WR&>(c);
}

#endif /* PARTICLE_CELL_H_ */
