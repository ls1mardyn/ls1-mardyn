/*
 * ParticleCellForwardDeclaration.h
 *
 *  Created on: 6 Feb 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_PARTICLECELLFORWARDDECLARATION_H_
#define SRC_PARTICLECONTAINER_PARTICLECELLFORWARDDECLARATION_H_

#ifndef ENABLE_REDUCED_MEMORY_MODE
	class FullParticleCell;
	typedef FullParticleCell ParticleCell;
#else
	class ParticleCellRMM;
	typedef ParticleCellRMM ParticleCell;
#endif

#endif /* SRC_PARTICLECONTAINER_PARTICLECELLFORWARDDECLARATION_H_ */
