#ifndef PARTICLE_CELL_H_
#define PARTICLE_CELL_H_

/**
 * the old class ParticleCell is now called FullParticleCell and it implements ParticleCellBase.
 */

#include "FullParticleCell.h"
#include "ParticleCellRMM.h"

#ifndef ENABLE_REDUCED_MEMORY_MODE
	typedef FullParticleCell ParticleCell;
#else
	typedef ParticleCellRMM ParticleCell;
#endif

inline FullParticleCell* downcastCellPointerFull(ParticleCellBase* c) {
	mardyn_assert(static_cast<FullParticleCell*>(c) == dynamic_cast<FullParticleCell*>(c));
	return static_cast<FullParticleCell*>(c);
}

inline FullParticleCell& downcastCellReferenceFull(ParticleCellBase& c) {
	mardyn_assert(&static_cast<FullParticleCell&>(c) == &dynamic_cast<FullParticleCell&>(c));
	return static_cast<FullParticleCell&>(c);
}

inline ParticleCellRMM* downcastCellPointerRMM(ParticleCellBase* c) {
	mardyn_assert(static_cast<ParticleCellRMM*>(c) == dynamic_cast<ParticleCellRMM*>(c));
	return static_cast<ParticleCellRMM*>(c);
}

inline ParticleCellRMM& downcastCellReferenceRMM(ParticleCellBase& c) {
	mardyn_assert(&static_cast<ParticleCellRMM&>(c) == &dynamic_cast<ParticleCellRMM&>(c));
	return static_cast<ParticleCellRMM&>(c);
}

#endif /* PARTICLE_CELL_H_ */
