#include "Molecule.h"
#include "particleContainer/adapter/CellDataSoA.h"

#include <cassert>
#include <cmath>
#include <fstream>

#include "utils/Logger.h"

using namespace std;
using Log::global_log;


Molecule::Molecule(unsigned long id, Component *component,
	                 double rx,  double ry,  double rz,
	                 double vx,  double vy,  double vz,
	                 double q0,  double q1,  double q2, double q3,
	                 double Dx,  double Dy,  double Dz
		  )
		: _q(q0, q1, q2, q3) {
	_id = id;
	_component = component;
	_r[0] = rx;
	_r[1] = ry;
	_r[2] = rz;
	_v[0] = vx;
	_v[1] = vy;
	_v[2] = vz;
	_L[0] = Dx;
	_L[1] = Dy;
	_L[2] = Dz;
	_Vi[0]= 0.;
	_Vi[1]= 0.;
	_Vi[2]= 0.;

	_soa = nullptr;
	_soa_index_lj = 0;
	_soa_index_c = 0;
	_soa_index_d = 0;
	_soa_index_q = 0;

	if(_component != nullptr) {
		_m = _component->m();
		_I[0] = _component->I11();
		_I[1] = _component->I22();
		_I[2] = _component->I33();
		for (unsigned short d = 0; d < 3; ++d) {
			if (_I[d] != 0.)
				_invI[d] = 1. / _I[d];
			else
				_invI[d] = 0.;
		}
	}
	_F[0] = _F[1] = _F[2] = 0.;
	_M[0] = _M[1] = _M[2] = 0.;
}

Molecule::Molecule(const Molecule& m) {
	_id = m._id;
	_component = m._component;
	_r[0] = m._r[0];
	_r[1] = m._r[1];
	_r[2] = m._r[2];
	_v[0] = m._v[0];
	_v[1] = m._v[1];
	_v[2] = m._v[2];
	_q = m._q;
	_L[0] = m._L[0];
	_L[1] = m._L[1];
	_L[2] = m._L[2];
	_F[0] = m._F[0];
	_F[1] = m._F[1];
	_F[2] = m._F[2];
	_M[0] = m._M[0];
	_M[1] = m._M[1];
	_M[2] = m._M[2];
	_Vi[0]= m._Vi[0];
	_Vi[1]= m._Vi[1];
	_Vi[2]= m._Vi[2];

	_soa = m._soa;
	_soa_index_lj = m._soa_index_lj;
	_soa_index_c = m._soa_index_c;
	_soa_index_d = m._soa_index_d;
	_soa_index_q = m._soa_index_q;

	_m = m._m;
	_I[0] = m._I[0];
	_I[1] = m._I[1];
	_I[2] = m._I[2];
	_invI[0] = m._invI[0];
	_invI[1] = m._invI[1];
	_invI[2] = m._invI[2];
}

std::array<double, 3> Molecule::ljcenter_d_abs(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_lj;
	std::array<double,3> ret;
	double* rx = _soa->ljc_r_xBegin();
	double* ry = _soa->ljc_r_yBegin();
	double* rz = _soa->ljc_r_zBegin();
	ret[0] = rx[index_in_soa];
	ret[1] = ry[index_in_soa];
	ret[2] = rz[index_in_soa];
	return ret;
}

std::array<double, 3> Molecule::charge_d_abs(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_c;
	std::array<double, 3> ret;
	double* rx = _soa->charges_r_xBegin();
	double* ry = _soa->charges_r_yBegin();
	double* rz = _soa->charges_r_zBegin();
	ret[0] = rx[index_in_soa];
	ret[1] = ry[index_in_soa];
	ret[2] = rz[index_in_soa];
	return ret;
}

std::array<double, 3> Molecule::dipole_d_abs(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_d;
	std::array<double, 3> ret;
	double* rx = _soa->dipoles_r_xBegin();
	double* ry = _soa->dipoles_r_yBegin();
	double* rz = _soa->dipoles_r_zBegin();
	ret[0] = rx[index_in_soa];
	ret[1] = ry[index_in_soa];
	ret[2] = rz[index_in_soa];
	return ret;
}

std::array<double, 3> Molecule::quadrupole_d_abs(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_q;
	std::array<double, 3> ret;
	double* rx = _soa->quadrupoles_r_xBegin();
	double* ry = _soa->quadrupoles_r_yBegin();
	double* rz = _soa->quadrupoles_r_zBegin();
	ret[0] = rx[index_in_soa];
	ret[1] = ry[index_in_soa];
	ret[2] = rz[index_in_soa];
	return ret;
}

std::array<double, 3> Molecule::ljcenter_F(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_lj;
	std::array<double, 3> ret;
	double* fx = _soa->ljc_f_xBegin();
	double* fy = _soa->ljc_f_yBegin();
	double* fz = _soa->ljc_f_zBegin();
	ret[0] = fx[index_in_soa];
	ret[1] = fy[index_in_soa];
	ret[2] = fz[index_in_soa];
	return ret;
}
std::array<double, 3> Molecule::charge_F(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_c;
	std::array<double, 3> ret;
	double* fx = _soa->charges_f_xBegin();
	double* fy = _soa->charges_f_yBegin();
	double* fz = _soa->charges_f_zBegin();
	ret[0] = fx[index_in_soa];
	ret[1] = fy[index_in_soa];
	ret[2] = fz[index_in_soa];
	return ret;
}
std::array<double, 3> Molecule::dipole_F(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_d;
	std::array<double, 3> ret;
	double* fx = _soa->dipoles_f_xBegin();
	double* fy = _soa->dipoles_f_yBegin();
	double* fz = _soa->dipoles_f_zBegin();
	ret[0] = fx[index_in_soa];
	ret[1] = fy[index_in_soa];
	ret[2] = fz[index_in_soa];
	return ret;
}
std::array<double, 3> Molecule::quadrupole_F(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_q;
	std::array<double, 3> ret;
	double* fx = _soa->quadrupoles_f_xBegin();
	double* fy = _soa->quadrupoles_f_yBegin();
	double* fz = _soa->quadrupoles_f_zBegin();
	ret[0] = fx[index_in_soa];
	ret[1] = fy[index_in_soa];
	ret[2] = fz[index_in_soa];
	return ret;
}

std::array<double, 3> Molecule::dipole_e(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_d;
	std::array<double, 3> ret;
	ret[0] = _soa->_dipoles_e.x(index_in_soa);
	ret[1] = _soa->_dipoles_e.y(index_in_soa);
	ret[2] = _soa->_dipoles_e.z(index_in_soa);
	return ret;
}
std::array<double, 3> Molecule::quadrupole_e(unsigned int i) const {
	const unsigned index_in_soa = i + _soa_index_q;
	std::array<double, 3> ret;
	ret[0] = _soa->_quadrupoles_e.x(index_in_soa);
	ret[1] = _soa->_quadrupoles_e.y(index_in_soa);
	ret[2] = _soa->_quadrupoles_e.z(index_in_soa);
	return ret;
}

void Molecule::Fljcenteradd(unsigned int i, double a[]) {
	const unsigned index_in_soa = i + _soa_index_lj;
	double* fx = _soa->ljc_f_xBegin();
	double* fy = _soa->ljc_f_yBegin();
	double* fz = _soa->ljc_f_zBegin();
	fx[index_in_soa] += a[0];
	fy[index_in_soa] += a[1];
	fz[index_in_soa] += a[2];
}

void Molecule::Fchargeadd(unsigned int i, double a[]) {
	const unsigned index_in_soa = i + _soa_index_c;
	double* fx = _soa->charges_f_xBegin();
	double* fy = _soa->charges_f_yBegin();
	double* fz = _soa->charges_f_zBegin();
	fx[index_in_soa] += a[0];
	fy[index_in_soa] += a[1];
	fz[index_in_soa] += a[2];
}

void Molecule::Fdipoleadd(unsigned int i, double a[]) {
	const unsigned index_in_soa = i + _soa_index_d;
	double* fx = _soa->dipoles_f_xBegin();
	double* fy = _soa->dipoles_f_yBegin();
	double* fz = _soa->dipoles_f_zBegin();
	fx[index_in_soa] += a[0];
	fy[index_in_soa] += a[1];
	fz[index_in_soa] += a[2];
}

void Molecule::Fquadrupoleadd(unsigned int i, double a[]) {
	const unsigned index_in_soa = i + _soa_index_q;
	double* fx = _soa->quadrupoles_f_xBegin();
	double* fy = _soa->quadrupoles_f_yBegin();
	double* fz = _soa->quadrupoles_f_zBegin();
	fx[index_in_soa] += a[0];
	fy[index_in_soa] += a[1];
	fz[index_in_soa] += a[2];
}

void Molecule::Fljcentersub(unsigned int i, double a[]) {
	double minusA[3] = {-a[0], -a[1], -a[2]};
	Fljcenteradd(i, minusA);
}
void Molecule::Fchargesub(unsigned int i, double a[]) {
	double minusA[3] = {-a[0], -a[1], -a[2]};
	Fchargeadd(i, minusA);
}
void Molecule::Fdipolesub(unsigned int i, double a[]) {
	double minusA[3] = {-a[0], -a[1], -a[2]};
	Fdipoleadd(i, minusA);
}
void Molecule::Fquadrupolesub(unsigned int i, double a[]) {
	double minusA[3] = {-a[0], -a[1], -a[2]};
	Fquadrupoleadd(i, minusA);
}

void Molecule::upd_preF(double dt) {
	assert(_m > 0);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / _m;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _F[d];
		_r[d] += dt * _v[d];
	}

	double w[3];
	_q.rotateinv(_L, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qhalfstep;
	_q.differentiate(w, qhalfstep);
	qhalfstep.scale(dt_halve);
	qhalfstep.add(_q);
	double qcorr = 1. / sqrt(qhalfstep.magnitude2());
	qhalfstep.scale(qcorr);
	for (unsigned short d = 0; d < 3; ++d)
		_L[d] += dt_halve * _M[d];
	qhalfstep.rotateinv(_L, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qincr;
	qhalfstep.differentiate(w, qincr);
	qincr.scale(dt);
	_q.add(qincr);
	qcorr = 1. / sqrt(_q.magnitude2());
	_q.scale(qcorr);

}

void Molecule::upd_postF(double dt_halve, double& summv2, double& sumIw2) {
	using std::isnan; // C++11 needed

	calcFM();

	double dtInv2m = dt_halve / _m;
	double v2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _F[d];
		v2 += _v[d] * _v[d];
		_L[d] += dt_halve * _M[d];
	}
    assert(!isnan(v2)); // catches NaN
    summv2 += _m * v2;

	double w[3];
	_q.rotateinv(_L, w); // L = D = Iw
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
    assert(!isnan(Iw2)); // catches NaN
	sumIw2 += Iw2;
}


double Molecule::U_rot() {
	double w[3];
	_q.rotateinv(_L, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	return 0.5 * Iw2;
}

void Molecule::calculate_mv2_Iw2(double& summv2, double& sumIw2) {
	summv2 += _m * v2();
	double w[3];
	_q.rotateinv(_L, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void Molecule::calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) {
	double vcx = _v[0] - offx;
	double vcy = _v[1] - offy;
	double vcz = _v[2] - offz;
	summv2 += _m * (vcx*vcx + vcy*vcy + vcz*vcz);

	double w[3];
	_q.rotateinv(_L, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void Molecule::scale_v(double s, double offx, double offy, double offz) {
	this->vsub(offx, offy, offz);
	this->scale_v(s);
	this->vadd(offx, offy, offz);
}

void Molecule::write(ostream& ostrm) const {
	ostrm << _id << "\t" << (_component->ID() + 1) << "\t"
	      << _r[0] << " " << _r[1] << " " << _r[2] << "\t"
	      << _v[0] << " " << _v[1] << " " << _v[2] << "\t"
	      << _q.qw() << " " << _q.qx() << " " << _q.qy() << " " << _q.qz() << "\t"
	      << _L[0] << " " << _L[1] << " " << _L[2] << "\t"
	      << _Vi[0] << " " << _Vi[1] << " " << _Vi[2] << "\t"
	      << endl;
}

// private functions
// these are only used when compiling molecule.cpp and therefore might be inlined without any problems

void Molecule::clearFM() {
	assert(_soa != nullptr);
	_F[0] = _F[1] = _F[2] = 0.;
	_M[0] = _M[1] = _M[2] = 0.;
	_Vi[0]= _Vi[1]= _Vi[2]= 0.;

	// clear SoA-cache (quickest way)
	unsigned ns = numLJcenters();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_lj;
		_soa->ljc_f_xBegin()[index_in_soa] = 0.0;
		_soa->ljc_f_yBegin()[index_in_soa] = 0.0;
		_soa->ljc_f_zBegin()[index_in_soa] = 0.0;
		_soa->ljc_V_xBegin()[index_in_soa] = 0.0;
		_soa->ljc_V_yBegin()[index_in_soa] = 0.0;
		_soa->ljc_V_zBegin()[index_in_soa] = 0.0;
	}
	ns = numCharges();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_c;
		_soa->charges_f_xBegin()[index_in_soa] = 0.0;
		_soa->charges_f_yBegin()[index_in_soa] = 0.0;
		_soa->charges_f_zBegin()[index_in_soa] = 0.0;
		_soa->charges_V_xBegin()[index_in_soa] = 0.0;
		_soa->charges_V_yBegin()[index_in_soa] = 0.0;
		_soa->charges_V_zBegin()[index_in_soa] = 0.0;
	}
	ns = numDipoles();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_d;
		_soa->dipoles_f_xBegin()[index_in_soa] = 0.0;
		_soa->dipoles_f_yBegin()[index_in_soa] = 0.0;
		_soa->dipoles_f_zBegin()[index_in_soa] = 0.0;
		_soa->dipoles_V_xBegin()[index_in_soa] = 0.0;
		_soa->dipoles_V_yBegin()[index_in_soa] = 0.0;
		_soa->dipoles_V_zBegin()[index_in_soa] = 0.0;
		_soa->_dipoles_M.x(index_in_soa) = 0.0;
		_soa->_dipoles_M.y(index_in_soa) = 0.0;
		_soa->_dipoles_M.z(index_in_soa) = 0.0;
	}
	ns = numQuadrupoles();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_q;
		_soa->quadrupoles_f_xBegin()[index_in_soa] = 0.0;
		_soa->quadrupoles_f_yBegin()[index_in_soa] = 0.0;
		_soa->quadrupoles_f_zBegin()[index_in_soa] = 0.0;
		_soa->quadrupoles_V_xBegin()[index_in_soa] = 0.0;
		_soa->quadrupoles_V_yBegin()[index_in_soa] = 0.0;
		_soa->quadrupoles_V_zBegin()[index_in_soa] = 0.0;
		_soa->_quadrupoles_M.x(index_in_soa) = 0.0;
		_soa->_quadrupoles_M.y(index_in_soa) = 0.0;
		_soa->_quadrupoles_M.z(index_in_soa) = 0.0;
	}
}

void Molecule::calcFM() {
	using std::isnan; // C++11 needed

	//_M[0] = _M[1] = _M[2] = 0.;
	unsigned int ns = numSites();
	for (unsigned int si = 0; si < ns; ++si) {
		const std::array<double,3> Fsite = site_F(si);
		const std::array<double,3> dsite = site_d(si);
#ifndef NDEBUG
		/*
		 * catches NaN assignments
		 */
		for (int d = 0; d < 3; d++) {
			if (isnan(dsite[d])) {
				global_log->error() << "Severe dsite[" << d << "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
			if (isnan(Fsite[d])) {
				global_log->error() << "Severe Fsite[" << d << "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
		}
#endif

		Fadd(Fsite.data());
		_M[0] += dsite[1] * Fsite[2] - dsite[2] * Fsite[1];
		_M[1] += dsite[2] * Fsite[0] - dsite[0] * Fsite[2];
		_M[2] += dsite[0] * Fsite[1] - dsite[1] * Fsite[0];
	}

	// accumulate virial, dipoles_M and quadrupoles_M:
	double temp_M[3] = { 0., 0., 0. };
	double temp_Vi[3] = { 0., 0., 0. };

	ns = numLJcenters();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_lj;
		temp_Vi[0] += _soa->ljc_V_xBegin()[index_in_soa];
		temp_Vi[1] += _soa->ljc_V_yBegin()[index_in_soa];
		temp_Vi[2] += _soa->ljc_V_zBegin()[index_in_soa];
	}
	ns = numCharges();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_c;
		temp_Vi[0] += _soa->charges_V_xBegin()[index_in_soa];
		temp_Vi[1] += _soa->charges_V_yBegin()[index_in_soa];
		temp_Vi[2] += _soa->charges_V_zBegin()[index_in_soa];
	}
	ns = numDipoles();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_d;
		temp_Vi[0] += _soa->dipoles_V_xBegin()[index_in_soa];
		temp_Vi[1] += _soa->dipoles_V_yBegin()[index_in_soa];
		temp_Vi[2] += _soa->dipoles_V_zBegin()[index_in_soa];
		temp_M[0] += _soa->_dipoles_M.x(index_in_soa);
		temp_M[1] += _soa->_dipoles_M.y(index_in_soa);
		temp_M[2] += _soa->_dipoles_M.z(index_in_soa);
	}
	ns = numQuadrupoles();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_q;
		temp_Vi[0] += _soa->quadrupoles_V_xBegin()[index_in_soa];
		temp_Vi[1] += _soa->quadrupoles_V_yBegin()[index_in_soa];
		temp_Vi[2] += _soa->quadrupoles_V_zBegin()[index_in_soa];
		temp_M[0] += _soa->_quadrupoles_M.x(index_in_soa);
		temp_M[1] += _soa->_quadrupoles_M.y(index_in_soa);
		temp_M[2] += _soa->_quadrupoles_M.z(index_in_soa);
	}
	temp_Vi[0] *= 0.5;
	temp_Vi[1] *= 0.5;
	temp_Vi[2] *= 0.5;
	assert(!isnan(temp_Vi[0]));
	assert(!isnan(temp_Vi[1]));
	assert(!isnan(temp_Vi[2]));
	Viadd(temp_Vi);
	Madd(temp_M);
}


/**
 * catches NaN values and missing data
 *
 * @note Use isnan from cmath to check for nan.
 * If that's not available (C99), compare the value with itself. If the value
 * is NaN, the comparison will evaluate to false (according to IEEE754 spec.)
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void Molecule::check(unsigned long id) {
#ifndef NDEBUG
	using std::isnan; // C++11 needed

	assert(_id == id);
	assert(_m > 0.0);
	for (int d = 0; d < 3; d++) {
		assert(!isnan(_r[d]));
		assert(!isnan(_v[d]));
		assert(!isnan(_L[d]));
		assert(!isnan(_F[d]));
		assert(!isnan(_M[d]));
		assert(!isnan(_I[d]));
		// assert(!isnan(_Vi[d]));
		assert(!isnan(_invI[d]));
	}
	if(isnan(_Vi[0]) || isnan(_Vi[1]) || isnan(_Vi[2]))
	{
	   cout << "\n\talert: molecule id " << id << " (internal cid " << this->_component->ID() << ") has virial _Vi = (" << _Vi[0] << ", " << _Vi[1] << ", " << _Vi[2] << ")";
	   _Vi[0] = 0.0;
	   _Vi[1] = 0.0;
	   _Vi[2] = 0.0;
	}
#endif
}
#pragma GCC diagnostic pop

bool Molecule::isLessThan(const Molecule& m2) const {
	if (_r[2] < m2.r(2))
		return true;
	else if (_r[2] > m2.r(2))
		return false;
	else {
		if (_r[1] < m2.r(1))
			return true;
		else if (_r[1] > m2.r(1))
			return false;
		else {
			if (_r[0] < m2.r(0))
				return true;
			else if (_r[0] > m2.r(0))
				return false;
			else {
				global_log->error() << "LinkedCells::isFirstParticle: both Particles have the same position" << endl;
				exit(1);
			}
		}
	}
	return false; /* Silence warnings about missing return statement */
}


std::ostream& operator<<( std::ostream& os, const Molecule& m ) {
	os << "ID: " << m.id() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	os << "F:  (" << m.F(0) << ", " << m.F(1) << ", " << m.F(2) << ")\n" ;
	os << "q:  [[" << m.q().qw() << ", " << m.q().qx() << "], [" << m.q().qy() << ", " << m.q().qz()<< "]]\n" ;
	os << "w:  (" << m.D(0) << ", " << m.D(1) << ", " << m.D(2) << ")\n";
	os << "Vi:  (" << m.Vi(0) << ", " << m.Vi(1) << ", " << m.Vi(2) << ")" ;
	return os;
}


unsigned long Molecule::totalMemsize() const {
	unsigned long size = sizeof (*this);
	return size;
}

void Molecule::setupSoACache(CellDataSoA* const s, unsigned iLJ, unsigned iC,
		unsigned iD, unsigned iQ) {
	setSoA(s);
	setStartIndexSoA_LJ(iLJ);
	setStartIndexSoA_C(iC);
	setStartIndexSoA_D(iD);
	setStartIndexSoA_Q(iQ);

	normalizeQuaternion();

	unsigned ns = numLJcenters();
	for (unsigned j = 0; j < ns; ++j) {
		double centerPos[3];
		computeLJcenter_d(j, centerPos);
		const unsigned ind = _soa_index_lj + j;
		_soa->ljc_m_r_xBegin()[ind] = _r[0];
		_soa->ljc_m_r_yBegin()[ind] = _r[1];
		_soa->ljc_m_r_zBegin()[ind] = _r[2];
		_soa->ljc_r_xBegin()[ind] = centerPos[0] + _r[0];
		_soa->ljc_r_yBegin()[ind] = centerPos[1] + _r[1];
		_soa->ljc_r_zBegin()[ind] = centerPos[2] + _r[2];
		_soa->_ljc_id[ind] = getComponentLookUpID() + j;
	}
	ns = numCharges();
	for (unsigned j = 0; j < ns; ++j) {
		double centerPos[3];
		computeCharge_d(j, centerPos);
		const unsigned ind = _soa_index_c + j;
		_soa->charges_m_r_xBegin()[ind] = _r[0];
		_soa->charges_m_r_yBegin()[ind] = _r[1];
		_soa->charges_m_r_zBegin()[ind] = _r[2];
		_soa->charges_r_xBegin()[ind] = centerPos[0] + _r[0];
		_soa->charges_r_yBegin()[ind] = centerPos[1] + _r[1];
		_soa->charges_r_zBegin()[ind] = centerPos[2] + _r[2];
		_soa->_charges_q[ind] = component()->charge(j).q();
	}
	ns = numDipoles();
	for (unsigned j = 0; j < ns; ++j) {
		double centerPos[3];
		computeDipole_d(j, centerPos);
		double orientation[3];
		computeDipole_e(j, orientation);
		const unsigned ind = _soa_index_d + j;
		_soa->dipoles_m_r_xBegin()[ind] = _r[0];
		_soa->dipoles_m_r_yBegin()[ind] = _r[1];
		_soa->dipoles_m_r_zBegin()[ind] = _r[2];
		_soa->dipoles_r_xBegin()[ind] = centerPos[0] + _r[0];
		_soa->dipoles_r_yBegin()[ind] = centerPos[1] + _r[1];
		_soa->dipoles_r_zBegin()[ind] = centerPos[2] + _r[2];
		_soa->_dipoles_p[ind] = component()->dipole(j).absMy();
		_soa->_dipoles_e.x(ind) = orientation[0];
		_soa->_dipoles_e.y(ind) = orientation[1];
		_soa->_dipoles_e.z(ind) = orientation[2];
	}
	ns = numQuadrupoles();
	for (unsigned j = 0; j < ns; ++j) {
		double centerPos[3];
		computeQuadrupole_d(j, centerPos);
		double orientation[3];
		computeQuadrupole_e(j, orientation);
		const unsigned ind = _soa_index_q + j;
		_soa->quadrupoles_m_r_xBegin()[ind] = _r[0];
		_soa->quadrupoles_m_r_yBegin()[ind] = _r[1];
		_soa->quadrupoles_m_r_zBegin()[ind] = _r[2];
		_soa->quadrupoles_r_xBegin()[ind] = centerPos[0] + _r[0];
		_soa->quadrupoles_r_yBegin()[ind] = centerPos[1] + _r[1];
		_soa->quadrupoles_r_zBegin()[ind] = centerPos[2] + _r[2];
		_soa->_quadrupoles_m[ind] = component()->quadrupole(j).absQ();
		_soa->_quadrupoles_e.x(ind) = orientation[0];
		_soa->_quadrupoles_e.y(ind) = orientation[1];
		_soa->_quadrupoles_e.z(ind) = orientation[2];
	}
}
