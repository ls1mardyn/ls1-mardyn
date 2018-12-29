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

	CanonicalEnsemble() :
			_N(0), _V(0), _T(0), _mu(0), _p(0), _E(0), _E_trans(0), _E_rot(0) {
	    _type = "NVT";
	}

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

	void updateGlobalVariable(ParticleContainer *particleContainer, GlobalVariable variable) override;




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
