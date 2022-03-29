#include "FullMolecule.h"
#include "particleContainer/adapter/CellDataSoA.h"

#include "utils/mardyn_assert.h"
#include <cmath>
#include <fstream>

#include "utils/Logger.h"

using namespace std;
using Log::global_log;


FullMolecule::FullMolecule(unsigned long id, Component *component,
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
	_Vi[3]= 0.;
	_Vi[4]= 0.;
	_Vi[5]= 0.;
	_Vi[6]= 0.;
	_Vi[7]= 0.;
	_Vi[8]= 0.;
	_upot = 0;
	_ViConstCorr = 0;
	_upotConstCorr = 0;

	_soa = nullptr;
	_soa_index_lj = 0;
	_soa_index_c = 0;
	_soa_index_d = 0;
	_soa_index_q = 0;

	this->updateMassInertia();

	_F[0] = _F[1] = _F[2] = 0.;
	_M[0] = _M[1] = _M[2] = 0.;
}

FullMolecule::FullMolecule(const FullMolecule& m) {
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
	_Vi[3]= m._Vi[3];
	_Vi[4]= m._Vi[4];
	_Vi[5]= m._Vi[5];
	_Vi[6]= m._Vi[6];
	_Vi[7]= m._Vi[7];
	_Vi[8]= m._Vi[8];
	_upot = m._upot;
	_ViConstCorr = m._ViConstCorr;
	_upotConstCorr = m._upotConstCorr;

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

FullMolecule& FullMolecule::operator=(const FullMolecule& m) {
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
	_Vi[3]= m._Vi[3];
	_Vi[4]= m._Vi[4];
	_Vi[5]= m._Vi[5];
	_Vi[6]= m._Vi[6];
	_Vi[7]= m._Vi[7];
	_Vi[8]= m._Vi[8];
	_upot = m._upot;
	_ViConstCorr = m._ViConstCorr;
	_upotConstCorr = m._upotConstCorr;

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

	return *this;
}

std::array<double, 3> FullMolecule::ljcenter_d_abs(unsigned int i) const {
	mardyn_assert(i < numLJcenters());

	const unsigned index_in_soa = i + _soa_index_lj;

	std::array<double,3> ret;
	std::array<vcp_real_calc, 3> temp;

	temp = _soa->getTripletCalc(CellDataSoA::QuantityType::CENTER_POSITION, ConcSites::SiteType::LJC, index_in_soa);
	ret[0] = static_cast<double>(temp[0]);
	ret[1] = static_cast<double>(temp[1]);
	ret[2] = static_cast<double>(temp[2]);

	return ret;
}

std::array<double, 3> FullMolecule::charge_d_abs(unsigned int i) const {
	mardyn_assert(i < numCharges());

	const unsigned index_in_soa = i + _soa_index_c;
	std::array<double, 3> ret;
	std::array<vcp_real_calc, 3> temp;

	temp = _soa->getTripletCalc(CellDataSoA::QuantityType::CENTER_POSITION, ConcSites::SiteType::CHARGE, index_in_soa);
	ret[0] = static_cast<double>(temp[0]);
	ret[1] = static_cast<double>(temp[1]);
	ret[2] = static_cast<double>(temp[2]);

	return ret;
}

std::array<double, 3> FullMolecule::dipole_d_abs(unsigned int i) const {
	mardyn_assert(i < numDipoles());

	const unsigned index_in_soa = i + _soa_index_d;
	std::array<double, 3> ret;
	std::array<vcp_real_calc, 3> temp;

	temp = _soa->getTripletCalc(CellDataSoA::QuantityType::CENTER_POSITION, ConcSites::SiteType::DIPOLE, index_in_soa);
	ret[0] = static_cast<double>(temp[0]);
	ret[1] = static_cast<double>(temp[1]);
	ret[2] = static_cast<double>(temp[2]);

	return ret;
}

std::array<double, 3> FullMolecule::quadrupole_d_abs(unsigned int i) const {
	mardyn_assert(i < numQuadrupoles());

	const unsigned index_in_soa = i + _soa_index_q;
	std::array<double, 3> ret;
	std::array<vcp_real_calc, 3> temp;

	temp = _soa->getTripletCalc(CellDataSoA::QuantityType::CENTER_POSITION, ConcSites::SiteType::QUADRUPOLE, index_in_soa);
	ret[0] = static_cast<double>(temp[0]);
	ret[1] = static_cast<double>(temp[1]);
	ret[2] = static_cast<double>(temp[2]);

	return ret;
}

std::array<double, 3> FullMolecule::ljcenter_F(unsigned int i) const {
	mardyn_assert(i < numLJcenters());

	const unsigned index_in_soa = i + _soa_index_lj;
	std::array<double, 3> ret;
	std::array<vcp_real_accum, 3> temp;

	temp = _soa->getTripletAccum(CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::LJC, index_in_soa);
	ret[0] = static_cast<double>(temp[0]);
	ret[1] = static_cast<double>(temp[1]);
	ret[2] = static_cast<double>(temp[2]);

	return ret;
}
std::array<double, 3> FullMolecule::charge_F(unsigned int i) const {
	mardyn_assert(i < numCharges());

	const unsigned index_in_soa = i + _soa_index_c;
	std::array<double, 3> ret;
	std::array<vcp_real_accum, 3> temp;

	temp = _soa->getTripletAccum(CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::CHARGE, index_in_soa);
	ret[0] = static_cast<double>(temp[0]);
	ret[1] = static_cast<double>(temp[1]);
	ret[2] = static_cast<double>(temp[2]);

	return ret;
}
std::array<double, 3> FullMolecule::dipole_F(unsigned int i) const {
	mardyn_assert(i < numDipoles());

	const unsigned index_in_soa = i + _soa_index_d;
	std::array<double, 3> ret;
	std::array<vcp_real_accum, 3> temp;

	temp = _soa->getTripletAccum(CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::DIPOLE, index_in_soa);
	ret[0] = static_cast<double>(temp[0]);
	ret[1] = static_cast<double>(temp[1]);
	ret[2] = static_cast<double>(temp[2]);

	return ret;
}
std::array<double, 3> FullMolecule::quadrupole_F(unsigned int i) const {
	mardyn_assert(i < numQuadrupoles());

	const unsigned index_in_soa = i + _soa_index_q;
	std::array<double, 3> ret;
	std::array<vcp_real_accum, 3> temp;

	temp = _soa->getTripletAccum(CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::QUADRUPOLE, index_in_soa);
	ret[0] = static_cast<double>(temp[0]);
	ret[1] = static_cast<double>(temp[1]);
	ret[2] = static_cast<double>(temp[2]);

	return ret;
}

std::array<double, 3> FullMolecule::dipole_e(unsigned int i) const {
	mardyn_assert(i < numDipoles());

	const unsigned index_in_soa = i + _soa_index_d;
	std::array<double, 3> ret;
	ret[0] = static_cast<double>(_soa->_dipoles_e.x(index_in_soa));
	ret[1] = static_cast<double>(_soa->_dipoles_e.y(index_in_soa));
	ret[2] = static_cast<double>(_soa->_dipoles_e.z(index_in_soa));
	return ret;
}
std::array<double, 3> FullMolecule::quadrupole_e(unsigned int i) const {
	mardyn_assert(i < numQuadrupoles());

	const unsigned index_in_soa = i + _soa_index_q;
	std::array<double, 3> ret;
	ret[0] = static_cast<double>(_soa->_quadrupoles_e.x(index_in_soa));
	ret[1] = static_cast<double>(_soa->_quadrupoles_e.y(index_in_soa));
	ret[2] = static_cast<double>(_soa->_quadrupoles_e.z(index_in_soa));
	return ret;
}

void FullMolecule::Fljcenteradd(unsigned int i, double a[]) {
	mardyn_assert(i < numLJcenters());

	const unsigned index_in_soa = i + _soa_index_lj;
	std::array<vcp_real_accum, 3> temp;

	temp = _soa->getTripletAccum(CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::LJC, index_in_soa);
	temp[0] += a[0];
	temp[1] += a[1];
	temp[2] += a[2];
	_soa->setTripletAccum(temp, CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::LJC, index_in_soa);
}

void FullMolecule::Fchargeadd(unsigned int i, double a[]) {
	mardyn_assert(i < numCharges());

	const unsigned index_in_soa = i + _soa_index_c;
	std::array<vcp_real_accum, 3> temp;

	temp = _soa->getTripletAccum(CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::CHARGE, index_in_soa);
	temp[0] += a[0];
	temp[1] += a[1];
	temp[2] += a[2];
	_soa->setTripletAccum(temp, CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::CHARGE, index_in_soa);
}

void FullMolecule::Fdipoleadd(unsigned int i, double a[]) {
	mardyn_assert(i < numDipoles());

	const unsigned index_in_soa = i + _soa_index_d;
	std::array<vcp_real_accum, 3> temp;

	temp = _soa->getTripletAccum(CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::DIPOLE, index_in_soa);
	temp[0] += a[0];
	temp[1] += a[1];
	temp[2] += a[2];
	_soa->setTripletAccum(temp, CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::DIPOLE, index_in_soa);
}

void FullMolecule::Fquadrupoleadd(unsigned int i, double a[]) {
	mardyn_assert(i < numQuadrupoles());

	const unsigned index_in_soa = i + _soa_index_q;
	std::array<vcp_real_accum, 3> temp;

	temp = _soa->getTripletAccum(CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::QUADRUPOLE, index_in_soa);
	temp[0] += a[0];
	temp[1] += a[1];
	temp[2] += a[2];
	_soa->setTripletAccum(temp, CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::QUADRUPOLE, index_in_soa);
}

void FullMolecule::Fljcentersub(unsigned int i, double a[]) {
	double minusA[3] = {-a[0], -a[1], -a[2]};
	Fljcenteradd(i, minusA);
}
void FullMolecule::Fchargesub(unsigned int i, double a[]) {
	double minusA[3] = {-a[0], -a[1], -a[2]};
	Fchargeadd(i, minusA);
}
void FullMolecule::Fdipolesub(unsigned int i, double a[]) {
	double minusA[3] = {-a[0], -a[1], -a[2]};
	Fdipoleadd(i, minusA);
}
void FullMolecule::Fquadrupolesub(unsigned int i, double a[]) {
	double minusA[3] = {-a[0], -a[1], -a[2]};
	Fquadrupoleadd(i, minusA);
}

void FullMolecule::upd_preF(double dt) {
	mardyn_assert(_m > 0);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / _m;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _F[d];
		_r[d] += dt * _v[d];
	}

	std::array<double, 3> w = _q.rotateinv(D_arr());
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
	w = qhalfstep.rotateinv(D_arr());
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qincr;
	qhalfstep.differentiate(w, qincr);
	qincr.scale(dt);
	_q.add(qincr);
	qcorr = 1. / sqrt(_q.magnitude2());
	_q.scale(qcorr);

}

void FullMolecule::upd_postF(double dt_halve, double& summv2, double& sumIw2) {
	using std::isnan; // C++11 needed

	//calcFM(); //NOTE: This was moved to simulation.cpp calculateForces() and is called in simulate()

	double dtInv2m = dt_halve / _m;
	double v2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		_v[d] += dtInv2m * _F[d];
		v2 += _v[d] * _v[d];
		_L[d] += dt_halve * _M[d];
	}
    mardyn_assert(!isnan(v2)); // catches NaN
    summv2 += _m * v2;

	std::array<double, 3> w = _q.rotateinv(D_arr()); // L = D = Iw
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
    mardyn_assert(!isnan(Iw2)); // catches NaN
	sumIw2 += Iw2;
}


double FullMolecule::U_rot() {
	std::array<double, 3> w = _q.rotateinv(D_arr());
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	return 0.5 * Iw2;
}

double FullMolecule::U_rot_2() {
	std::array<double, 3> w = _q.rotateinv(D_arr());
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	return Iw2;
}

void FullMolecule::calculate_mv2_Iw2(double& summv2, double& sumIw2) {
	summv2 += _m * v2();
	std::array<double, 3> w = _q.rotateinv(D_arr());
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void FullMolecule::calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) {
	double vcx = _v[0] - offx;
	double vcy = _v[1] - offy;
	double vcz = _v[2] - offz;
	summv2 += _m * (vcx*vcx + vcy*vcy + vcz*vcz);

	std::array<double, 3> w = _q.rotateinv(D_arr());
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

std::string FullMolecule::getWriteFormat(){
	return std::string("ICRVQD");
}

void FullMolecule::write(ostream& ostrm) const {
	ostrm << _id << "\t" << (_component->ID() + 1) << "\t"
	      << _r[0] << " " << _r[1] << " " << _r[2] << "\t"
	      << _v[0] << " " << _v[1] << " " << _v[2] << "\t"
	      << _q.qw() << " " << _q.qx() << " " << _q.qy() << " " << _q.qz() << "\t"
	      << _L[0] << " " << _L[1] << " " << _L[2] << "\t"
	      << "\n";
}

void FullMolecule::writeBinary(std::ostream& ostrm) const {
	unsigned int cid = _component->ID() + 1;
	double qw = _q.qw();
	double qx = _q.qx();
	double qy = _q.qy();
	double qz = _q.qz();

	ostrm.write(reinterpret_cast<const char*>(&_id), 8);
	ostrm.write(reinterpret_cast<const char*>(&cid), 4);
	ostrm.write(reinterpret_cast<const char*>(&(_r[0])), 8);
	ostrm.write(reinterpret_cast<const char*>(&(_r[1])), 8);
	ostrm.write(reinterpret_cast<const char*>(&(_r[2])), 8);
	ostrm.write(reinterpret_cast<const char*>(&(_v[0])), 8);
	ostrm.write(reinterpret_cast<const char*>(&(_v[1])), 8);
	ostrm.write(reinterpret_cast<const char*>(&(_v[2])), 8);
	ostrm.write(reinterpret_cast<const char*>(&qw), 8);
	ostrm.write(reinterpret_cast<const char*>(&qx), 8);
	ostrm.write(reinterpret_cast<const char*>(&qy), 8);
	ostrm.write(reinterpret_cast<const char*>(&qz), 8);
	ostrm.write(reinterpret_cast<const char*>(&(_L[0])), 8);
	ostrm.write(reinterpret_cast<const char*>(&(_L[1])), 8);
	ostrm.write(reinterpret_cast<const char*>(&(_L[2])), 8);
}


// private functions
// these are only used when compiling molecule.cpp and therefore might be inlined without any problems

void FullMolecule::clearFM() {
	mardyn_assert(_soa != nullptr);
	_F[0] = _F[1] = _F[2] = 0.;
	_M[0] = _M[1] = _M[2] = 0.;
	_Vi[0]= _Vi[1]= _Vi[2]=_Vi[3]= _Vi[4]= _Vi[5]=_Vi[6]= _Vi[7]= _Vi[8]= 0.;
	_upot = 0.;

	std::array<vcp_real_accum, 3> clearance = {0.0, 0.0, 0.0};

	// clear SoA-cache (quickest way)
	unsigned ns = numLJcenters();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_lj;

		_soa->setTripletAccum(clearance, CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::LJC, index_in_soa);
		_soa->setTripletAccum(clearance, CellDataSoA::QuantityType::VIRIAL, ConcSites::SiteType::LJC, index_in_soa);
	}
	ns = numCharges();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_c;

		_soa->setTripletAccum(clearance, CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::CHARGE, index_in_soa);
		_soa->setTripletAccum(clearance, CellDataSoA::QuantityType::VIRIAL, ConcSites::SiteType::CHARGE, index_in_soa);
	}
	ns = numDipoles();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_d;

		_soa->setTripletAccum(clearance, CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::DIPOLE, index_in_soa);
		_soa->setTripletAccum(clearance, CellDataSoA::QuantityType::VIRIAL, ConcSites::SiteType::DIPOLE, index_in_soa);

		_soa->_dipoles_M.x(index_in_soa) = 0.0;
		_soa->_dipoles_M.y(index_in_soa) = 0.0;
		_soa->_dipoles_M.z(index_in_soa) = 0.0;
	}
	ns = numQuadrupoles();
	for (unsigned i = 0; i < ns; ++i) {
		const unsigned index_in_soa = i + _soa_index_q;

		_soa->setTripletAccum(clearance, CellDataSoA::QuantityType::FORCE, ConcSites::SiteType::QUADRUPOLE, index_in_soa);
		_soa->setTripletAccum(clearance, CellDataSoA::QuantityType::VIRIAL, ConcSites::SiteType::QUADRUPOLE, index_in_soa);

		_soa->_quadrupoles_M.x(index_in_soa) = 0.0;
		_soa->_quadrupoles_M.y(index_in_soa) = 0.0;
		_soa->_quadrupoles_M.z(index_in_soa) = 0.0;
	}
}

void FullMolecule::calcFM_site(const std::array<double, 3>& dsite, const std::array<double, 3>& Fsite) {
#ifndef NDEBUG
	using std::isnan; // C++11 needed

	/*
	 * catches NaN assignments
	 */
	for (int d = 0; d < 3; d++) {
		if (isnan(dsite[d])) {
			global_log->error() << "Severe dsite[" << d << "] error for site of m" << _id << endl;
			mardyn_assert(false);
		}
		if (isnan(Fsite[d])) {
			global_log->error() << "Severe Fsite[" << d << "] error for site of m" << _id << endl;
			mardyn_assert(false);
		}
	}
#endif

	Fadd(Fsite.data());
	_M[0] += dsite[1] * Fsite[2] - dsite[2] * Fsite[1];
	_M[1] += dsite[2] * Fsite[0] - dsite[0] * Fsite[2];
	_M[2] += dsite[0] * Fsite[1] - dsite[1] * Fsite[0];
}

void FullMolecule::calcFM() {
	using std::isnan; // C++11 needed

	//_M[0] = _M[1] = _M[2] = 0.;
	unsigned int ns;

	// accumulate virial, dipoles_M and quadrupoles_M:
	double temp_M[3] = { 0., 0., 0. };
	double temp_Vi[3] = { 0., 0., 0. };

	std::array<vcp_real_accum, 3> interim;

	ns = numLJcenters();
	for (unsigned i = 0; i < ns; ++i) {
		const std::array<double,3> Fsite = ljcenter_F(i);
		const std::array<double,3> dsite = ljcenter_d(i);
		calcFM_site(dsite, Fsite);

		const unsigned index_in_soa = i + _soa_index_lj;
		interim = _soa->getTripletAccum(CellDataSoA::QuantityType::VIRIAL, ConcSites::SiteType::LJC, index_in_soa);

		temp_Vi[0] += interim[0];
		temp_Vi[1] += interim[1];
		temp_Vi[2] += interim[2];
	}
	ns = numCharges();
	for (unsigned i = 0; i < ns; ++i) {
		const std::array<double,3> Fsite = charge_F(i);
		const std::array<double,3> dsite = charge_d(i);
		calcFM_site(dsite, Fsite);

		const unsigned index_in_soa = i + _soa_index_c;
		interim = _soa->getTripletAccum(CellDataSoA::QuantityType::VIRIAL, ConcSites::SiteType::CHARGE, index_in_soa);

		temp_Vi[0] += interim[0];
		temp_Vi[1] += interim[1];
		temp_Vi[2] += interim[2];
	}
	ns = numDipoles();
	for (unsigned i = 0; i < ns; ++i) {
		const std::array<double,3> Fsite = dipole_F(i);
		const std::array<double,3> dsite = dipole_d(i);
		calcFM_site(dsite, Fsite);

		const unsigned index_in_soa = i + _soa_index_d;
		interim = _soa->getTripletAccum(CellDataSoA::QuantityType::VIRIAL, ConcSites::SiteType::DIPOLE, index_in_soa);

		temp_Vi[0] += interim[0];
		temp_Vi[1] += interim[1];
		temp_Vi[2] += interim[2];
		temp_M[0] += _soa->_dipoles_M.x(index_in_soa);
		temp_M[1] += _soa->_dipoles_M.y(index_in_soa);
		temp_M[2] += _soa->_dipoles_M.z(index_in_soa);
	}
	ns = numQuadrupoles();
	for (unsigned i = 0; i < ns; ++i) {
		const std::array<double,3> Fsite = quadrupole_F(i);
		const std::array<double,3> dsite = quadrupole_d(i);
		calcFM_site(dsite, Fsite);

		const unsigned index_in_soa = i + _soa_index_q;
		interim = _soa->getTripletAccum(CellDataSoA::QuantityType::VIRIAL, ConcSites::SiteType::QUADRUPOLE, index_in_soa);

		temp_Vi[0] += interim[0];
		temp_Vi[1] += interim[1];
		temp_Vi[2] += interim[2];
		temp_M[0] += _soa->_quadrupoles_M.x(index_in_soa);
		temp_M[1] += _soa->_quadrupoles_M.y(index_in_soa);
		temp_M[2] += _soa->_quadrupoles_M.z(index_in_soa);
	}
	temp_Vi[0] *= 0.5;
	temp_Vi[1] *= 0.5;
	temp_Vi[2] *= 0.5;
	mardyn_assert(!isnan(temp_Vi[0]));
	mardyn_assert(!isnan(temp_Vi[1]));
	mardyn_assert(!isnan(temp_Vi[2]));
	Viadd(temp_Vi);
	Madd(temp_M);
}


/**
 * catches NaN values and missing data
 */
#ifndef NDEBUG
void FullMolecule::check(unsigned long id) {

  using std::isfinite; // C++11 needed

  mardyn_assert(_id == id);
  mardyn_assert(_m > 0.0);
  for (int d = 0; d < 3; d++) {
    mardyn_assert(isfinite(_r[d]));
    mardyn_assert(isfinite(_v[d]));
    mardyn_assert(isfinite(_L[d]));
    mardyn_assert(isfinite(_F[d]));
    mardyn_assert(isfinite(_M[d]));
    mardyn_assert(isfinite(_I[d]));
    // mardyn_assert(!isnan(_Vi[d]));
    mardyn_assert(isfinite(_invI[d]));
  }
  _q.check();
  if (!isfinite(_Vi[0]) || !isfinite(_Vi[1]) || !isfinite(_Vi[2]
   || !isfinite(_Vi[3]) || !isfinite(_Vi[4]) || !isfinite(_Vi[5])
   || !isfinite(_Vi[6]) || !isfinite(_Vi[7]) || !isfinite(_Vi[8]))) {
    cout << "\talert: molecule id " << id << " (internal cid " << this->_component->ID() << ") has virial _Vi = ("
         << _Vi[0] << ", " << _Vi[1] << ", " << _Vi[2] << ", "
         << _Vi[3] << ", " << _Vi[4] << ", " << _Vi[5] << ", "
         << _Vi[6] << ", " << _Vi[7] << ", " << _Vi[8] << ")" << endl;
    _Vi[0] = 0.0;
    _Vi[1] = 0.0;
    _Vi[2] = 0.0;
    _Vi[3] = 0.0;
    _Vi[4] = 0.0;
    _Vi[5] = 0.0;
    _Vi[6] = 0.0;
    _Vi[7] = 0.0;
    _Vi[8] = 0.0;
    mardyn_assert(false);
  }
}
#else
void FullMolecule::check(unsigned long /*id*/) {}
#endif

std::ostream& operator<<( std::ostream& os, const FullMolecule& m ) {
	os << "ID: " << m.getID() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	os << "F:  (" << m.F(0) << ", " << m.F(1) << ", " << m.F(2) << ")\n" ;
	os << "q:  [[" << m.q().qw() << ", " << m.q().qx() << "], [" << m.q().qy() << ", " << m.q().qz()<< "]]\n" ;
	os << "w:  (" << m.D(0) << ", " << m.D(1) << ", " << m.D(2) << ")\n";
	os << "Vi:  (" << m.Vi(0) << ", " << m.Vi(3) << ", " << m.Vi(4) << ")" ;  // Elements written out as in tensor
	os << "Vi:  (" << m.Vi(6) << ", " << m.Vi(1) << ", " << m.Vi(5) << ")" ;  // Elements written out as in tensor
	os << "Vi:  (" << m.Vi(7) << ", " << m.Vi(8) << ", " << m.Vi(2) << ")" ;  // Elements written out as in tensor
	return os;
}


unsigned long FullMolecule::totalMemsize() const {
	unsigned long size = sizeof (*this);
	return size;
}

void FullMolecule::setSoA(CellDataSoABase * const s) {
	CellDataSoA * derived;
	derived = static_cast<CellDataSoA *>(s);
	_soa = derived;
}

void FullMolecule::buildOwnSoA() {
	const size_t nLJ = numLJcenters();
	const size_t nC = numCharges();
	const size_t nD = numDipoles();
	const size_t nQ = numQuadrupoles();

	_soa = new CellDataSoA(1ul, nLJ, nC, nD, nQ);

	_soa->_mol_ljc_num[0] = nLJ;
	_soa->_mol_charges_num[0] = nC;
	_soa->_mol_dipoles_num[0] = nD;
	_soa->_mol_quadrupoles_num[0] = nQ;

	_soa->_mol_pos.x(0) = r(0);
	_soa->_mol_pos.y(0) = r(1);
	_soa->_mol_pos.z(0) = r(2);

	setupSoACache(_soa, 0ul, 0ul, 0ul, 0ul);
	clearFM();
}

void FullMolecule::releaseOwnSoA() {
	delete _soa;
}

void FullMolecule::setupSoACache(CellDataSoABase* const s, unsigned iLJ, unsigned iC,
		unsigned iD, unsigned iQ) {
	setSoA(s);
	setStartIndexSoA_LJ(iLJ);
	setStartIndexSoA_C(iC);
	setStartIndexSoA_D(iD);
	setStartIndexSoA_Q(iQ);

	normalizeQuaternion();

	unsigned ns = numLJcenters();
	for (unsigned j = 0; j < ns; ++j) {
		std::array<double, 3> centerPos = computeLJcenter_d(j);
		centerPos[0] += _r[0];
		centerPos[1] += _r[1];
		centerPos[2] += _r[2];

		const unsigned ind = _soa_index_lj + j;

		_soa->pushBackLJC(ind, convert_double_to_vcp_real_calc(r_arr()), convert_double_to_vcp_real_calc(centerPos), getComponentLookUpID() + j);
	}
	ns = numCharges();
	for (unsigned j = 0; j < ns; ++j) {
		std::array<double, 3> centerPos = computeCharge_d(j);
		centerPos[0] += _r[0];
		centerPos[1] += _r[1];
		centerPos[2] += _r[2];

		const unsigned ind = _soa_index_c + j;

		_soa->pushBackCharge(ind, convert_double_to_vcp_real_calc(r_arr()), convert_double_to_vcp_real_calc(centerPos), component()->charge(j).q());
	}
	ns = numDipoles();
	for (unsigned j = 0; j < ns; ++j) {
		std::array<double, 3> centerPos = computeDipole_d(j);
		centerPos[0] += _r[0];
		centerPos[1] += _r[1];
		centerPos[2] += _r[2];

		std::array<double,3> orientation = computeDipole_e(j);
		const unsigned ind = _soa_index_d + j;

		_soa->pushBackDipole(ind, convert_double_to_vcp_real_calc(r_arr()), convert_double_to_vcp_real_calc(centerPos), component()->dipole(j).absMy(), convert_double_to_vcp_real_calc(orientation));
	}
	ns = numQuadrupoles();
	for (unsigned j = 0; j < ns; ++j) {
		std::array<double, 3> centerPos = computeQuadrupole_d(j);
		centerPos[0] += _r[0];
		centerPos[1] += _r[1];
		centerPos[2] += _r[2];

		std::array<double,3> orientation = computeQuadrupole_e(j);
		const unsigned ind = _soa_index_q + j;

		_soa->pushBackQuadrupole(ind, convert_double_to_vcp_real_calc(r_arr()), convert_double_to_vcp_real_calc(centerPos), component()->quadrupole(j).absQ(), convert_double_to_vcp_real_calc(orientation));
	}
}
