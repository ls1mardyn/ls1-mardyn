#ifndef CANONICAL_ENSEMBLE_H_
#define CANONICAL_ENSEMBLE_H_

#include <vector>

#include "ensemble/EnsembleBase.h"
#include "particleContainer/ParticleContainer.h"

/** @brief Canonical ensemble (NVT) 
 *  @author Christoph Niethammer <niethammer@hlrs.de>
 *  
 * This class provides access to all global variables of the canonical ensemble (NVT).
 **/
class Component;

class DomainBase;
/* Fix problem with Cray compiler which requires the actual size of the component class. */
#include "molecules/Component.h"

class CanonicalEnsemble : public Ensemble {

private:
	CanonicalEnsemble& operator=(CanonicalEnsemble ensemble);

public:

	CanonicalEnsemble();

	virtual ~CanonicalEnsemble() {
	}

	virtual void readXML(XMLfileUnits& xmlconfig) override;

	unsigned long N() override {
		return _N;
	}

	double V() override {
		return _V;
	}

	double T() override {
		return _T;
	}

	double mu() override {
		return _mu;
	}

	double p() override {
		return _p;
	}

	double E() override {
		return _E;
	}

	void updateGlobalVariable(std::vector<ParticleContainer*>& particleContainers, GlobalVariable variable) override;

	/*! runs before temperature control is applied, but after force calculations, triggers radial distribution call in Domain only for NVT */
	void beforeThermostat(unsigned long simstep, unsigned long initStatistics) override;

private:

	int numComponents() {
		return _components.size();
	}

	unsigned long _N;
	double _V;
	double _T;

	double _mu;
	double _p;
	double _E;

	double _E_trans;
	double _E_rot;

	// TODO: As the canonical ensemble fixes the temperature here should be the right place for the thermostats
	//std::map<int, VelocityScalingThermostat> _thermostats;
	//std::map<int, int> _componentToThermostatMap;

};

#endif  /* CANONICAL_ENSEMBLE_H_ */
