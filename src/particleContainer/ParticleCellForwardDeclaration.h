/*
 * ParticleCellForwardDeclaration.h
 *
 *  Created on: 6 Feb 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_PARTICLECELLFORWARDDECLARATION_H_
#define SRC_PARTICLECONTAINER_PARTICLECELLFORWARDDECLARATION_H_

#ifndef MARDYN_WR
	class FullParticleCell;
	typedef FullParticleCell ParticleCell;
#else
	class ParticleCell_WR;
	typedef ParticleCell_WR ParticleCell;
#endif

#endif /* SRC_PARTICLECONTAINER_PARTICLECELLFORWARDDECLARATION_H_ */
