/*
 * MultipoleParticle.h
 *
 *  Created on: Nov 27, 2014
 *      Author: tchipevn
 */

#ifndef MULTIPOLEPARTICLE_H_
#define MULTIPOLEPARTICLE_H_

#include "bhfmm/pseudoParticles/PseudoParticle.h"

namespace bhfmm {
class LocalParticle;

/**
 * Interface for Multipole particles of any expansion.
 */
class MultipoleParticle: public PseudoParticle {
public:
	MultipoleParticle() :
			PseudoParticle() {
	}

	virtual ~MultipoleParticle() {
	}

	/**
	 * P2M operator
	 * @param position
	 * @param charge
	 */
	virtual void addSource(const Vector3<double>& position, double charge) = 0;

	/**
	 * M2M operator
	 * @param small the smaller particle to be added to the larger one (this)
	 */
	virtual void addMultipoleParticle(const MultipoleParticle& small) = 0;

	/**
	 * M2L operator
	 * @param local
	 */
	virtual void actOnLocalParticle(LocalParticle& local) const = 0;

	/**
	 * M2P operator
	 * @param position
	 * @param charge
	 * @param potential stores resulting potential
	 * @param force stores resulting force
	 */
	virtual void actOnTarget(const Vector3<double>& position, double charge, double& pot,
			Vector3<double>& force) const = 0;

};

} /* namespace bhfmm */

#endif /* MULTIPOLEPARTICLE_H_ */
