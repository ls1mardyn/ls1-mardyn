/***************************************************************************
 *   Copyright (C) 2009 by Martin Bernreuther and colleagues               *
 *   Martin.Bernreuther@informatik.uni-stuttgart.de                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 *   Due to copyleft all future versions of this program must be           *
 *   distributed as Free Software (e.g., using a BSD-like license).        *
 ***************************************************************************/
#include <cmath>
#include <fstream>
#include <cassert>

#include "Molecule.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

Molecule::Molecule(unsigned long id, unsigned int componentid, double rx,
		double ry, double rz, double vx, double vy, double vz, double q0,
		double q1, double q2, double q3, double Dx, double Dy, double Dz,
		const vector<Component>* components) :
	_q(q0, q1, q2, q3), _ljcenters(NULL), _tersoff(NULL) {
	_id = id;
	_componentid = componentid;
	_r[0] = rx;
	_r[1] = ry;
	_r[2] = rz;
	_v[0] = vx;
	_v[1] = vy;
	_v[2] = vz;
	_D[0] = Dx;
	_D[1] = Dy;
	_D[2] = Dz;
	_sites_d = _sites_F = _osites_e = NULL;
	_numTersoffNeighbours = 0;
	fixedx = rx;
	fixedy = ry;

	_leftxSite = new double[6];
	_leftxF = new double[3];
	_leftxRdfSite = new double[6];
	_leftxRdfF = new double[3];

	for (int i = 0; i < 3; i++) {
		_leftxF[i] = 0;
		_leftxRdfF[i] = 0;
	}
	for (int i = 0; i < 6; i++) {
		_leftxRdfSite[i] = 0;
		_leftxSite[i] = 0;
	}

	if (components)
		setupCache(components);

	_low_boundary = _high_boundary = -1;
	_low_boundary_speed = _high_boundary_speed = 0;
}

Molecule::Molecule(const Molecule& m) {

	_leftxSite = new double[6];
	_leftxF = new double[3];
	_leftxRdfSite = new double[6];
	_leftxRdfF = new double[3];

	for (int i = 0; i < 3; i++) {
		_leftxRdfF[i] = 0;
		_leftxF[i] = 0;
	}
	for (int i = 0; i < 6; i++) {
		_leftxRdfSite[i] = 0;
		_leftxSite[i] = 0;
	}

	_id = m._id;
	_componentid = m._componentid;
	_r[0] = m._r[0];
	_r[1] = m._r[1];
	_r[2] = m._r[2];
	_v[0] = m._v[0];
	_v[1] = m._v[1];
	_v[2] = m._v[2];
	_q = m._q;
	_D[0] = m._D[0];
	_D[1] = m._D[1];
	_D[2] = m._D[2];
	_F[0] = m._F[0];
	_F[1] = m._F[1];
	_F[2] = m._F[2];
	_M[0] = m._M[0];
	_M[1] = m._M[1];
	_M[2] = m._M[2];

	_ljcenters = m._ljcenters;
	_charges = m._charges;
	_dipoles = m._dipoles;
	_quadrupoles = m._quadrupoles;
	assert( m._tersoff );
	_tersoff = m._tersoff;
	_m = m._m;
	_I[0] = m._I[0];
	_I[1] = m._I[1];
	_I[2] = m._I[2];
	_invI[0] = m._invI[0];
	_invI[1] = m._invI[1];
	_invI[2] = m._invI[2];

	_numsites = m._numsites;
	_numorientedsites = m._numorientedsites;
	assert(_numsites);
	_sites_d = new double[_numsites * 3];
	assert(_sites_d);

	for (unsigned int i = 0; i < _numsites * 3; ++i)
		_sites_d[i] = m._sites_d[i]; // not necessary -> cache only
	_ljcenters_d = &(_sites_d[0]);
	_charges_d = &(_ljcenters_d[numLJcenters() * 3]);
	_dipoles_d = &(_charges_d[numCharges() * 3]);
	_quadrupoles_d = &(_dipoles_d[numDipoles() * 3]);
	_tersoff_d = &(_quadrupoles_d[numQuadrupoles() * 3]);

	_osites_e = new double[_numorientedsites * 3];
	assert(_osites_e);
	//for(unsigned int i=0;i<_numorientedsites*3;++i) _osites_e[i]=m._osites_e[i]; // not necessary -> cache only
	_dipoles_e = &(_osites_e[0]);
	_quadrupoles_e = &(_dipoles_e[numDipoles() * 3]);

	_sites_F = new double[_numsites * 3];

	assert(_sites_F);
	//for(unsigned int i=0;i<_numsites*3;++i) _sites_F[i]=m._sites_F[i]; // not necessary -> cache only
	_ljcenters_F = &(_sites_F[0]);
	_charges_F = &(_ljcenters_F[numLJcenters() * 3]);
	_dipoles_F = &(_charges_F[numCharges() * 3]);
	_quadrupoles_F = &(_dipoles_F[numDipoles() * 3]);
	_tersoff_F = &(_quadrupoles_F[numQuadrupoles() * 3]);
	_numTersoffNeighbours = 0;
	fixedx = m.fixedx;
	fixedy = m.fixedy;

	_low_boundary = _high_boundary = -1;
	_low_boundary_speed = _high_boundary_speed = 0;
}

void Molecule::enableBouncingBack(double lowb, double highb, double lowspeed,
		double highspeed) {
	_low_boundary = lowb;
	_high_boundary = highb;
	_low_boundary_speed = lowspeed;
	_high_boundary_speed = highspeed;
}

void Molecule::upd_preF(double dt, double vcorr, double Dcorr) {
	assert(_m > 0);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / _m;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] = vcorr * _v[d] + dtInv2m * _F[d];
		_r[d] += dt * _v[d];

	}

	if ((_low_boundary != _high_boundary && _low_boundary != -1) && (_r[0]
			< _low_boundary || _r[0] > _high_boundary))
		bounceBackDirection(dt, 0);
	_dt = dt;

	double w[3];
	_q.rotate(_D, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qhalfstep;
	_q.differentiate(w, qhalfstep);
	qhalfstep.scale(dt_halve);
	qhalfstep.add(_q);
	double qcorr = 1. / sqrt(qhalfstep.magnitude2());
	if (isnan(qcorr)) {
		std::cout << "qcorr nan for " << id() << std::endl;
	}
	qhalfstep.scale(qcorr);
	for (unsigned short d = 0; d < 3; ++d)
		_D[d] = Dcorr * _D[d] + dt_halve * _M[d];
	qhalfstep.rotate(_D, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qincr;
	qhalfstep.differentiate(w, qincr);
	qincr.scale(dt);
	_q.add(qincr);
	qcorr = 1. / sqrt(_q.magnitude2());
	_q.scale(qcorr);

}

void Molecule::upd_cache() {
	unsigned int i;
	unsigned int ns;

	ns = numLJcenters();
	double mag = 1 / std::sqrt(_q.magnitude2());
	if (isnan(mag) && id() == 1) {
		std::cout << "velocity " << v(0) << " " << v(1) << " " << v(2)
				<< std::endl;
		std::cout << "angular velocity " << D(0) << " " << D(1) << " " << D(2)
				<< std::endl;
		std::cout << "force: " << F(0) << " " << F(1) << " " << F(2)
				<< std::endl;
		std::cout << "nan for molecule " << id() << std::endl;
	}
	_q.scale(mag);
	//_q *= std::sqrt(_q.magnitude2());
	for (i = 0; i < ns; ++i)
		_q.rotateinv((*_ljcenters)[i].r(), &(_ljcenters_d[i * 3]));
	ns = numCharges();
	for (i = 0; i < ns; ++i)
		_q.rotateinv((*_charges)[i].r(), &(_charges_d[i * 3]));
	ns = numDipoles();
	for (i = 0; i < ns; ++i) {
		const Dipole& di = (*_dipoles)[i];
		_q.rotateinv(di.r(), &(_dipoles_d[i * 3]));
		_q.rotateinv(di.e(), &(_dipoles_e[i * 3]));
	}
	ns = numQuadrupoles();
	for (i = 0; i < ns; ++i) {
		const Quadrupole& qi = (*_quadrupoles)[i];
		_q.rotateinv(qi.r(), &(_quadrupoles_d[i * 3]));
		_q.rotateinv(qi.e(), &(_quadrupoles_e[i * 3]));
	}
	ns = numTersoff();
	for (i = 0; i < ns; i++)
		_q.rotateinv((*_tersoff)[i].r(), &(_tersoff_d[i * 3]));
}

void Molecule::upd_postF(double dt_halve, double& summv2, double& sumIw2) {

	calcFM();

	double dtInv2m = dt_halve / _m;
	double v2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _F[d];
		v2 += _v[d] * _v[d];
		_D[d] += dt_halve * _M[d];
	}
	assert(!isnan(v2)); // catches NaN
	summv2 += _m * v2;

	double w[3];
	_q.rotate(_D, w); // L = D = Iw
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	assert(!isnan(Iw2)); // catches NaN
	sumIw2 += Iw2;

	/*
	 // if reflective wall should be used
	 if (_low_boundary == _high_boundary && _low_boundary == -1)
	 return;

	 // do reflective wall


	 // get the time of crossing the lower and upper boundary
	 double t_c_low = (_r[0] - _low_boundary) / (_low_boundary_speed - _v[0]);
	 double t_c_high = (_r[0] - _high_boundary) / (_high_boundary_speed - _v[0]);


	 // if crossing the lower boundary (2.1 for nummerical instabilities)
	 // boudnary speed is 0 for equilibrium state
	 if (t_c_low <= dt_halve * 2.1 && t_c_low > 0) {
	 double old_v = _v[0];

	 _v[0] = old_v - 2 * (old_v - _low_boundary_speed);
	 _r[0] += t_c_low * old_v + (2 * dt_halve - t_c_low) * _v[0];

	 // for trying random orientations
	 //		Quaternion q = Quaternion(rand(), rand(), rand(), rand());
	 //		q.normalize();
	 //
	 //		setq(q);
	 }

	 // if crossing higher boundary
	 if (t_c_high < dt_halve * 2.1 && t_c_high > 0) {
	 double old_v = _v[0];

	 _v[0] = old_v - 2 * (old_v - _high_boundary_speed);
	 _r[0] += t_c_high * old_v + (2 * dt_halve - t_c_high) * _v[0];


	 // for trying random orientations
	 //		Quaternion q = Quaternion(rand(), rand(), rand(), rand());
	 //		q.normalize();
	 //
	 //		setq(q);
	 }
	 */
}

double Molecule::U_rot() {
	double w[3];
	_q.rotate(_D, w);
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
	_q.rotate(_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void Molecule::calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx,
		double offy, double offz) {
	double vcx = _v[0] - offx;
	double vcy = _v[1] - offy;
	double vcz = _v[2] - offz;
	summv2 += _m * (vcx * vcx + vcy * vcy + vcz * vcz);

	double w[3];
	_q.rotate(_D, w);
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
	ostrm << _id << "\t" << (_componentid + 1) << "\t" << _r[0] << " " << _r[1]
			<< " " << _r[2] << "\t" << _v[0] << " " << _v[1] << " " << _v[2]
			<< "\t" << _q.qw() << " " << _q.qx() << " " << _q.qy() << " "
			<< _q.qz() << "\t" << _D[0] << " " << _D[1] << " " << _D[2] << "\t"
			<< endl;
}

void Molecule::addTersoffNeighbour(Molecule* m, bool pairType) {
	// this->_Tersoff_neighbours.insert(pair<Molecule*, bool>(m, (pairType > 0)));
	for (int j = 0; j < _numTersoffNeighbours; j++) {
		if (m->_id == _Tersoff_neighbours_first[j]->id()) {
			this->_Tersoff_neighbours_first[j] = m;
			this->_Tersoff_neighbours_second[j] = pairType;
			return;
		}
	}

	this->_Tersoff_neighbours_first[_numTersoffNeighbours] = m;
	this->_Tersoff_neighbours_second[_numTersoffNeighbours] = pairType;
	this->_numTersoffNeighbours++;
	if (_numTersoffNeighbours > MAX_TERSOFF_NEIGHBOURS) {
		global_log->error() << "Tersoff neighbour list overflow: Molecule "
				<< m->_id << " has more than " << MAX_TERSOFF_NEIGHBOURS
				<< " Tersoff neighbours." << endl;
		exit(1);
	}
}

double Molecule::tersoffParameters(double params[15]) //returns delta_r
{
	const Tersoff* t = &((*(this->_tersoff))[0]);
	params[0] = t->R();
	params[1] = t->S();
	params[2] = t->h();
	params[3] = t->cSquare();
	params[4] = t->dSquare();
	params[5] = t->A();
	params[6] = t->minusLambda();
	params[7] = t->minusMu();
	params[8] = t->beta();
	params[9] = t->n();
	params[10] = M_PI / (t->S() - t->R());
	params[11] = 1.0 + t->cSquare() / t->dSquare();
	params[12] = t->S() * t->S();
	params[13] = -(t->B());
	params[14] = -0.5 / t->n();

	return 0.000001 * (t->S() - t->R());
}

// private functions
// these are only used when compiling molecule.cpp and therefore might be inlined without any problems

inline void Molecule::setupCache(const vector<Component>* components) {
	assert(components);
	if (components->size() == 0)
		return;
	_numsites = _numorientedsites = 0;
	_ljcenters = &(*components)[_componentid].ljcenters();
	_numsites += _ljcenters->size();
	_charges = &(*components)[_componentid].charges();
	_numsites += _charges->size();
	_dipoles = &(*components)[_componentid].dipoles();
	_numsites += _dipoles->size();
	_numorientedsites += _dipoles->size();
	_quadrupoles = &(*components)[_componentid].quadrupoles();
	_numsites += _quadrupoles->size();
	_numorientedsites += _quadrupoles->size();
	_tersoff = &(*components)[_componentid].tersoff();
#ifndef NDEBUG
	if (!_tersoff) {
		global_log->error()
				<< "Tersoff vector null pointer detected for Molecule " << _id
				<< endl;
		exit(1);
	}
#endif
	_numsites += _tersoff->size();

	_m = (*components)[_componentid].m();
	_I[0] = (*components)[_componentid].I11();
	_I[1] = (*components)[_componentid].I22();
	_I[2] = (*components)[_componentid].I33();
	for (unsigned short d = 0; d < 3; ++d) {
		if (_I[d] != 0.)
			_invI[d] = 1. / _I[d];
		else
			_invI[d] = 0.;
	}

	assert(_numsites);

	_sites_d = new double[_numsites * 3];

	assert(_sites_d);
	_ljcenters_d = &(_sites_d[0]);
	_charges_d = &(_ljcenters_d[numLJcenters() * 3]);
	_dipoles_d = &(_charges_d[numCharges() * 3]);
	_quadrupoles_d = &(_dipoles_d[numDipoles() * 3]);
	_tersoff_d = &(_quadrupoles_d[numQuadrupoles() * 3]);

	_osites_e = new double[_numorientedsites * 3];
	assert(_osites_e);
	_dipoles_e = &(_osites_e[0]);
	_quadrupoles_e = &(_dipoles_e[numDipoles() * 3]);

	_sites_F = new double[_numsites * 3];

	assert(_sites_F);
	_ljcenters_F = &(_sites_F[0]);
	_charges_F = &(_ljcenters_F[numLJcenters() * 3]);
	_dipoles_F = &(_charges_F[numCharges() * 3]);
	_quadrupoles_F = &(_dipoles_F[numDipoles() * 3]);
	_tersoff_F = &(_quadrupoles_F[numQuadrupoles() * 3]);

	this->clearFM();
}

void Molecule::clearFM() {
	for (unsigned int i = 0; i < _numsites * 3; ++i)
		_sites_F[i] = 0.;
	_F[0] = _F[1] = _F[2] = 0.;
	_M[0] = _M[1] = _M[2] = 0.;
}

void Molecule::calcFM() {
	//_M[0] = _M[1] = _M[2] = 0.;
	unsigned int ns = numSites();
	for (unsigned int si = 0; si < ns; ++si) {
		const double* Fsite = site_F(si);
		const double* dsite = site_d(si);

#ifndef NDEBUG
		/*
		 * catches NaN assignments
		 */
		for (int d = 0; d < 3; d++) {
			if (isnan(dsite[d])) {
				global_log->error() << "Severe dsite[" << d
						<< "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
			if (isnan(Fsite[d])) {
				global_log->error() << "Severe Fsite[" << d
						<< "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
		}
#endif

		Fadd(Fsite);

		_M[0] += dsite[1] * Fsite[2] - dsite[2] * Fsite[1];
		_M[1] += dsite[2] * Fsite[0] - dsite[0] * Fsite[2];
		_M[2] += dsite[0] * Fsite[1] - dsite[1] * Fsite[0];
	}
}

// tijana
void Molecule::calcLeftxInfluence() {
	unsigned int ns = numSites();
	for (int i = 0; i < 3; i++) {
		_leftxF[i] = 0;
		_leftxRdfF[i] = 0;
	}
	for (unsigned int si = 0; si < ns; ++si) {
		const double* leftxSite = &_leftxSite[3 * si];
		const double* leftxRdfSite = &_leftxRdfSite[3 * si];
		leftxFAdd(leftxSite);
		leftxRdfFAdd(leftxRdfSite);
	}
}
/*
 * catches NaN values and missing data
 *
 * @note Use isnan from cmath to check for nan.
 * If that's not available (C99), compare the value with itself. If the value
 * is NaN, the comparison will evaluate to false (according to IEEE754 spec.)
 */
void Molecule::check(unsigned long id) {
#ifndef NDEBUG
	assert(_id == id);
	assert(_m > 0.0);
	assert(_numsites > 0);
	for (int d = 0; d < 3; d++) {
		assert(!isnan(_r[d]));
		assert(!isnan(_v[d]));
		assert(!isnan(_D[d]));
		assert(!isnan(_F[d]));
		assert(!isnan(_M[d]));
		assert(!isnan(_I[d]));
		assert(!isnan(_invI[d]));
	}
#endif
}

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
				global_log->error()
						<< "LinkedCells::isFirstParticle: both Particles have the same position"
						<< endl;
				exit(1);
			}
		}
	}
	return false; /* Silence warnings about missing return statement */
}

std::ostream& operator<<(std::ostream& os, const Molecule& m) {
	os << "ID: " << m.id() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n";
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n";
	os << "F:  (" << m.F(0) << ", " << m.F(1) << ", " << m.F(2) << ")\n";
	os << "q:  [[" << m.q().qw() << ", " << m.q().qx() << "], [" << m.q().qy()
			<< ", " << m.q().qz() << "]]\n";
	os << "w:  (" << m.D(0) << ", " << m.D(1) << ", " << m.D(2) << ")";
	return os;
}

unsigned long Molecule::totalMemsize() const {
	unsigned long size = sizeof(*this);

	//_sites_d
	size += sizeof(double) * _numsites * 3;
	// site orientation _osites_e
	size += sizeof(double) * _numorientedsites * 3;
	// site Forces _sites_F
	// row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
	size += sizeof(double) * _numsites * 3;

	return size;
}

void Molecule::setv(double* v) {
	for (int i = 0; i < 3; i++)
		_v[i] = v[i];
}

void Molecule::setD(double* D) {
	for (int i = 0; i < 3; i++) {
		_D[i] = D[i];
	}
}

void Molecule::bounceBack(double dt) {
	for (int d = 0; d < 3; d++) {
		_v[d] = -_v[d];
	}

	for (int d = 0; d < 3; d++) {
		_r[d] += dt * _v[d];

	}
}

void Molecule::bounceBackDirection(double dt, int dir) {
	_v[dir] = -_v[dir];
	_r[dir] += dt * _v[dir];
}
