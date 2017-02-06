#ifndef PARTICLE_CELL_H_
#define PARTICLE_CELL_H_

/**
 * the old class ParticleCell is now called FullParticleCell and it implements ParticleCellBase.
 */

#ifndef MARDYN_WR
	#include "FullParticleCell.h"
	typedef FullParticleCell ParticleCell;
#else
	#include "ParticleCellWR.h"
	typedef ParticleCell_WR ParticleCell;
#endif

#endif /* PARTICLE_CELL_H_ */
