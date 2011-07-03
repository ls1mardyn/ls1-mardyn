#ifndef CANONICAL_ENSEMBLE_H
#define CANONICAL_ENSEMBLE_H

#include <vector>
#include <map>

#include "ensemble/EnsembleBase.h"
#include "particleContainer/ParticleContainer.h"
//#include "ensemble/VelocityScalingThermostat.h"

/** @brief Canonical ensemble (NVT) 
 *  @author Christoph Niethammer <niethammer@hlrs.de>
 *  
 * This class provides access to all global variables of the canonical ensemble (NVT).
 **/
class Component;
/* Fix problem with Cray compiler which requires the actual size of the component class. */
#include "molecules/Component.h" 


class CanonicalEnsemble : public Ensemble {
public:

	CanonicalEnsemble() : _N(0), _V(0), _T(0), _mu(0), _p(0), _E(0), _E_trans(0), _E_rot(0),
		_components( 0 ), _particles( 0 )  {}

	CanonicalEnsemble( ParticleContainer *particles, std::vector<Component> *components) :
		_N(0), _V(0), _T(0), _mu(0), _p(0), _E(0), _E_trans(0), _E_rot(0),
		_components(components), _particles(particles)  {}
	~CanonicalEnsemble(){}

	unsigned long N() { return _N; }
	double V() { return _V; }
	double T() { return _T; }

	double mu(){ return _mu;}
	double p() { return _p; }
	double E() { return _E; }

	void updateGlobalVariable( GlobalVariable variable );


	int numComponents() { return _components->size(); }
	std::vector<Component>& components() { return *_components; }

private:
	unsigned long _N;
	double _V;
	double _T;

	double _mu;
	double _p;
	double _E;

	double _E_trans;
	double _E_rot;

	// cached pointers for easier access to particles and components,
	std::vector<Component> *_components;
	ParticleContainer *_particles;
	// As the canonical ensemble fixes the temperature here should be the right place for the thermostats
	//std::map<int, VelocityScalingThermostat> _thermostats;
	//std::map<int, int> _componentToThermostatMap;

};

#endif  // CANONICAL_ENSEMBLE_H
