/**
 * @file AutoPasFullMolecule.cpp
 * @author seckler
 * @date 20.09.18
 */

#include "AutoPasSimpleMolecule.h"
#include "Simulation.h"

Component* AutoPasSimpleMolecule::_component = nullptr;

Quaternion AutoPasSimpleMolecule::_quaternion = Quaternion(1.0, 0.0, 0.0, 0.0);

AutoPasSimpleMolecule::AutoPasSimpleMolecule(unsigned long id, Component* component, double rx, double ry, double rz,
											 double vx, double vy, double vz, double q0, double q1, double q2,
											 double q3, double Dx, double Dy, double Dz)
	: autopas::MoleculeLJ<double>({rx, ry, rz}, {vx, vy, vz}, id) {
	if (_component == nullptr) {
		_component = component;
	} else if (_component != component and component != nullptr) {
		global_log->debug() << "AutoPasSimpleMolecule can only handle one component" << std::endl;
		_component = component;
		// Simulation::exit(32);
	}
}

void AutoPasSimpleMolecule::upd_preF(double dt) {
	double mass = component()->m();
	mardyn_assert(mass > 0);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / mass;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _f[d];
		_r[d] += dt * _v[d];
	}
}

void AutoPasSimpleMolecule::upd_postF(double dt_halve, double& summv2, double& sumIw2) {
	using std::isnan;  // C++11 needed

	double mass = component()->m();
	double dtInv2m = dt_halve / mass;
	double v2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _f[d];
		v2 += _v[d] * _v[d];
	}
	mardyn_assert(!isnan(v2));  // catches NaN
	summv2 += mass * v2;
}
