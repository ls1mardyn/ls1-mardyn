/*
 * LocalParticle.h
 *
 *  Created on: Nov 27, 2014
 *      Author: tchipevn
 */

#ifndef LOCALPARTICLE_H_
#define LOCALPARTICLE_H_

#include "bhfmm/pseudoParticles/PseudoParticle.h"

namespace bhfmm {
class MultipoleParticle;

class LocalParticle: public PseudoParticle {
public:
	LocalParticle() :
			PseudoParticle() {
	}

	virtual ~LocalParticle() {
	}

	/**
	 * P2L operator
	 * @param position
	 * @param charge
	 */
	virtual void addSource(const Vector3<double>& position, double charge) = 0;

	/**
	 * M2L operator
	 * @param multipole
	 * @param periodicShift - extra shift for periodic boundary conditions
	 * @todo remove periodicShift, once the new container is up and running
	 */
	virtual void addMultipoleParticle(const MultipoleParticle& multipole, Vector3<double> periodicShift) = 0;

	/**
	 * L2L operator
	 * @param small
	 */
	virtual void actOnLocalParticle(LocalParticle& small) const = 0;

	/**
	 * L2P operator
	 * @param position
	 * @param charge
	 * @param potential stores resulting potential
	 * @param force stores resulting force
	 */
	virtual void actOnTarget(const Vector3<double>& position, double charge, double& potential,
			Vector3<double>& force) const = 0;
};

} /* namespace bhfmm */

#endif /* LOCALPARTICLE_H_ */
