/**
 * @file AutoPasFullMolecule.cpp
 * @author seckler
 * @date 20.09.18
 */

#include "AutoPasSimpleMolecule.h"
#include "utils/mardyn_assert.h"

Component* AutoPasSimpleMolecule::_component = nullptr;

Quaternion AutoPasSimpleMolecule::_quaternion = Quaternion(1.0, 0.0, 0.0, 0.0);

AutoPasSimpleMolecule::AutoPasSimpleMolecule(unsigned long id, Component* component, double rx, double ry, double rz,
											 double vx, double vy, double vz, double qw, double qx, double qy,
											 double qz, double Dx, double Dy, double Dz)
	: autopas::ParticleBaseFP64({rx, ry, rz}, {vx, vy, vz}, id) {
	if (_component == nullptr) {
		_component = component;
	} else if (_component != component and component != nullptr) {
		Log::global_log->warning() << "AutoPasSimpleMolecule can only handle one component" << std::endl;
		_component = component;
		// MARDYN_EXIT(error_message.str());
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
	mardyn_assert(!std::isnan(v2));  // catches NaN
	summv2 += mass * v2;
}

std::ostream& operator<<( std::ostream& os, const AutoPasSimpleMolecule& m ) {
	os << "ID: " << m.getID() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	os << "F:  (" << m.F(0) << ", " << m.F(1) << ", " << m.F(2) << ")\n" ;
	os << "q:  [[" << m.q().qw() << ", " << m.q().qx() << "], [" << m.q().qy() << ", " << m.q().qz()<< "]]\n" ;
	os << "w:  (" << m.D(0) << ", " << m.D(1) << ", " << m.D(2) << ")\n";
	os << "Vi:  (" << m.Vi(0) << ", " << m.Vi(1) << ", " << m.Vi(2) << ")" ;
	return os;
}
