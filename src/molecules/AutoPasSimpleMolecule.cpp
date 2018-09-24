/**
 * @file AutoPasFullMolecule.cpp
 * @author seckler
 * @date 20.09.18
 */

#include "AutoPasSimpleMolecule.h"

const std::array<double, 3> &AutoPasSimpleMolecule::getR() const {
	// this works as long as std::array is plain old data.
	auto &pos = reinterpret_cast<const std::array<double, 3> &>(_r);
	mardyn_assert(pos.data() == _r);
	return pos;
}

bool AutoPasSimpleMolecule::inBox(const std::array<double, 3> &rmin, const std::array<double, 3> &rmax) const {
	return MoleculeInterface::inBox(rmin.data(), rmax.data());
}

std::string AutoPasSimpleMolecule::toString() {
	std::ostringstream text;
	// clang-format off
	text << "Particle"
	     << "\nID      : " << _id
	     << "\nPosition: "
	     << _r[0] << " | " << _r[1] << " | " << _r[2]
	     << "\nVelocity: "
	     << _v[0] << " | " << _v[1] << " | " << _v[2]
	     << "\nForce   : "
	     << _F[0] << " | " << _F[1] << " | " << _F[2];
	// clang-format on
	return text.str();
}

const std::array<double, 3> &AutoPasSimpleMolecule::getF() const {
	// this works as long as std::array is plain old data.
	auto &force = reinterpret_cast<const std::array<double, 3> &>(_F);
	mardyn_assert(force.data() == _F);
	return force;
}

void AutoPasSimpleMolecule::setF(const std::array<double, 3>& F) {
	for (unsigned short i = 0; i < 3; i++) {
		_F[i] = F[i];
	}
}

void AutoPasSimpleMolecule::setF(double F[3]) {
	for (unsigned short i = 0; i < 3; i++) {
		_F[i] = F[i];
	}
}

void AutoPasSimpleMolecule::addF(const std::array<double, 3> &F) {
	for (unsigned short i = 0; i < 3; i++) {
		_F[i] += F[i];
	}
}

void AutoPasSimpleMolecule::subF(const std::array<double, 3> &F) {
	for (unsigned short i = 0; i < 3; i++) {
		_F[i] -= F[i];
	}
}