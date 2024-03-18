#ifndef SRC_INTEGRATORS_LEAPFROGRMM_H_
#define SRC_INTEGRATORS_LEAPFROGRMM_H_

#include "Integrator.h"

class VelocityCellProcessorRMM;

class LeapfrogRMM : public Integrator {
public:
	LeapfrogRMM();
	LeapfrogRMM (double timestepLength);
	~LeapfrogRMM();

	void readXML(XMLfileUnits& xmlconfig);

	void init() {}

	void eventForcesCalculated(ParticleContainer* moleculeContainer, Domain* domain) {
		computeVelocities(moleculeContainer, domain);
	}

	void eventNewTimestep(ParticleContainer* moleculeContainer, Domain* domain) {
		computePositions(moleculeContainer, domain);
	}

private:

	void computePositions(ParticleContainer* molCont, Domain* dom);
	void computeVelocities(ParticleContainer* molCont, Domain* dom);

	// unlike the PositionCellProcessor, the VelocityCellProcessor has dynamic data, so keep it.
	VelocityCellProcessorRMM * _velocityCellProcessor;

};

#endif /* SRC_INTEGRATORS_LEAPFROGRMM_H_ */

