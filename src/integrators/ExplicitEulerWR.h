/*
 * LeapFrogWR.h
 *
 *  Created on: Apr 16, 2017
 *      Author: tchipevn
 */

#ifndef SRC_INTEGRATORS_EXPLICITEULERWR_H_
#define SRC_INTEGRATORS_EXPLICITEULERWR_H_

#include "Integrator.h"

class ExplicitEuler_WR : public Integrator {
public:
	ExplicitEuler_WR() {}
	ExplicitEuler_WR(double timestepLength);
	~ExplicitEuler_WR() {}

	void readXML(XMLfileUnits& xmlconfig);

	void init() {}

	void eventForcesCalculated(ParticleContainer* moleculeContainer, Domain* domain) {
		computeVelocities(moleculeContainer, domain);
	}

	void eventNewTimestep(ParticleContainer* moleculeContainer, Domain* domain) {
		computePositions(moleculeContainer, domain);
	}

	void accelerateUniformly(
			ParticleContainer* molCont,
			Domain* domain
	) {}
	void accelerateInstantaneously(
			ParticleContainer* molCont,
			Domain* domain
	) {}

private:
	
	void computePositions(ParticleContainer* molCont, Domain* dom);
	void computeVelocities(ParticleContainer* molCont, Domain* dom);

};

#endif /* SRC_INTEGRATORS_EXPLICITEULERWR_H_ */
