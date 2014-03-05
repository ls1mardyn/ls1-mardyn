#include "Molecule.h"

#include <cassert>
#include <cmath>
#include <fstream>

#include "utils/Logger.h"

using namespace std;
using Log::global_log;


Molecule::Molecule(unsigned long id, Component *component,
	                 double rx, double ry, double rz,
	                 double vx, double vy, double vz,
	                 double q0, double q1, double q2, double q3,
	                 double Dx, double Dy, double Dz)
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
	_sites_d = _sites_F =_osites_e = NULL;
	_numTersoffNeighbours = 0;
	fixedx = rx;
	fixedy = ry;

	if(_component != NULL) {
		setupCache();
	}
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
	_sites_d = _sites_F =_osites_e = NULL;
	_numTersoffNeighbours = 0;
	fixedx = m.fixedx;
	fixedy = m.fixedy;

	if(_component != NULL) {
		setupCache();
	}
}


void Molecule::upd_preF(double dt, double vcorr, double Dcorr) {
	assert(_m > 0);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / _m;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] = vcorr * _v[d] + dtInv2m * _F[d];
		_r[d] += dt * _v[d];
	}

	double w[3];
	_q.rotate(_L, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qhalfstep;
	_q.differentiate(w, qhalfstep);
	qhalfstep.scale(dt_halve);
	qhalfstep.add(_q);
	double qcorr = 1. / sqrt(qhalfstep.magnitude2());
	qhalfstep.scale(qcorr);
	for (unsigned short d = 0; d < 3; ++d)
		_L[d] = Dcorr * _L[d] + dt_halve * _M[d];
	qhalfstep.rotate(_L, w);
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

	_q.normalize();

	ns = numLJcenters();
	for (i = 0; i < ns; ++i)
		_q.rotateinv(_component->ljcenter(i).r(), &(_ljcenters_d[i*3]));
	ns = numCharges();
	for (i = 0; i < ns; ++i)
		_q.rotateinv(_component->charge(i).r(), &(_charges_d[i*3]));
	ns = numDipoles();
	for (i = 0; i < ns; ++i) {
		const Dipole& di = _component->dipole(i);
		_q.rotateinv(di.r(), &(_dipoles_d[i*3]));
		_q.rotateinv(di.e(), &(_dipoles_e[i*3]));
	}
	ns = numQuadrupoles();
	for (i = 0; i < ns; ++i) {
		const Quadrupole& qi = _component->quadrupole(i);
		_q.rotateinv(qi.r(), &(_quadrupoles_d[i*3]));
		_q.rotateinv(qi.e(), &(_quadrupoles_e[i*3]));
	}
	ns = numTersoff();
	for (i = 0; i < ns; i++)
		_q.rotateinv(_component->tersoff(i).r(), &(_tersoff_d[i*3]));
}


void Molecule::upd_postF(double dt_halve, double& summv2, double& sumIw2) {
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
	_q.rotate(_L, w); // L = D = Iw
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
	_q.rotate(_L, w);
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
	_q.rotate(_L, w);
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
	_q.rotate(_L, w);
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
		global_log->error() << "Tersoff neighbour list overflow: Molecule " << m->_id << " has more than " << MAX_TERSOFF_NEIGHBOURS << " Tersoff neighbours." << endl;
		exit(1);
	}
}

double Molecule::tersoffParameters(double params[15]) //returns delta_r
{
	const Tersoff* t = &_component->tersoff()[0];
	params[ 0] = t->R();
	params[ 1] = t->S();
	params[ 2] = t->h();
	params[ 3] = t->cSquare();
	params[ 4] = t->dSquare();
	params[ 5] = t->A();
	params[ 6] = t->minusLambda();
	params[ 7] = t->minusMu();
	params[ 8] = t->beta();
	params[ 9] = t->n();
	params[10] = M_PI / (t->S() - t->R());
	params[11] = 1.0 + t->cSquare() / t->dSquare();
	params[12] = t->S() * t->S();
	params[13] = -(t->B());
	params[14] = -0.5 / t->n();

	return 0.000001 * (t->S() - t->R());
}

// private functions
// these are only used when compiling molecule.cpp and therefore might be inlined without any problems

inline void Molecule::setupCache() {
	assert(_component != NULL);

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

	int numsites = _component->numSites();
	_sites_d = new double[3*numsites];
	assert(_sites_d);
	_ljcenters_d = &(_sites_d[0]);
	_charges_d = &(_ljcenters_d[3*numLJcenters()]);
	_dipoles_d = &(_charges_d[3*numCharges()]);
	_quadrupoles_d = &(_dipoles_d[3*numDipoles()]);
	_tersoff_d = &(_quadrupoles_d[3*numQuadrupoles()]);

	int numorientedsites = _component->numOrientedSites();
	_osites_e = new double[3*numorientedsites];
	assert(_osites_e);
	_dipoles_e = &(_osites_e[0]);
	_quadrupoles_e = &(_dipoles_e[3*numDipoles()]);

	_sites_F = new double[3*numsites];
	assert(_sites_F);
	_ljcenters_F = &(_sites_F[0]);
	_charges_F = &(_ljcenters_F[3*numLJcenters()]);
	_dipoles_F = &(_charges_F[3*numCharges()]);
	_quadrupoles_F = &(_dipoles_F[3*numDipoles()]);
	_tersoff_F = &(_quadrupoles_F[3*numQuadrupoles()]);

	this->clearFM();
}

void Molecule::clearFM() {
	int numSites = _component->numSites();
	for (int i = 0; i < 3*numSites; i++) {
		_sites_F[i] = 0.;
	}
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
				global_log->error() << "Severe dsite[" << d << "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
			if (isnan(Fsite[d])) {
				global_log->error() << "Severe Fsite[" << d << "] error for site " << si << " of m" << _id << endl;
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


/**
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
	for (int d = 0; d < 3; d++) {
		assert(!isnan(_r[d]));
		assert(!isnan(_v[d]));
		assert(!isnan(_L[d]));
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
	os << "w:  (" << m.D(0) << ", " << m.D(1) << ", " << m.D(2) << ")" ;
	return os;
}


unsigned long Molecule::totalMemsize() const {
	unsigned long size = sizeof (*this);

	//_sites_d
	size += sizeof(double) * _component->numSites()* 3;
	// site orientation _osites_e
	size += sizeof(double) * _component->numOrientedSites() * 4;
	// site Forces _sites_F
	// row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
	size += sizeof(double) * _component->numSites() * 3;

	return size;
}
