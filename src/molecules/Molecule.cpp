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
#include "molecules/Molecule.h"


#include <cmath>
#include <fstream>
#include <cassert>

#include "utils/Logger.h"

using namespace std;
using Log::global_log;

Domain* Molecule::_domain;

Molecule::Molecule(unsigned long id, int componentid,
	                 double rx, double ry, double rz,
	                 double vx, double vy, double vz,
	                 double q0, double q1, double q2, double q3,
	                 double Dx, double Dy, double Dz,
	                 const vector<Component>* components)
		: m_q(q0, q1, q2, q3), m_ljcenters(NULL), m_tersoff(NULL) {
	m_id = id;
	m_componentid = componentid;
	m_r[0] = rx;
	m_r[1] = ry;
	m_r[2] = rz;
	m_v[0] = vx;
	m_v[1] = vy;
	m_v[2] = vz;
	m_D[0] = Dx;
	m_D[1] = Dy;
	m_D[2] = Dz;
	m_sites_d = m_sites_F = m_osites_e = NULL;
	m_curTN = 0;
	fixedx = rx;
	fixedy = ry;
	if (components)
		setupCache(components);
}

Molecule::Molecule(const Molecule& m) {
	m_id = m.m_id;
	m_componentid = m.m_componentid;
	m_r[0] = m.m_r[0];
	m_r[1] = m.m_r[1];
	m_r[2] = m.m_r[2];
	m_v[0] = m.m_v[0];
	m_v[1] = m.m_v[1];
	m_v[2] = m.m_v[2];
	m_q = m.m_q;
	m_D[0] = m.m_D[0];
	m_D[1] = m.m_D[1];
	m_D[2] = m.m_D[2];
	m_F[0] = m.m_F[0];
	m_F[1] = m.m_F[1];
	m_F[2] = m.m_F[2];
	m_M[0] = m.m_M[0];
	m_M[1] = m.m_M[1];
	m_M[2] = m.m_M[2];

	m_ljcenters = m.m_ljcenters;
	m_charges = m.m_charges;
	m_dipoles = m.m_dipoles;
	m_quadrupoles = m.m_quadrupoles;
	if (!m.m_tersoff) {
		global_log->warning() << "Tersoff vector null pointer detected for Molecule " << m_id << endl;
	}
	m_tersoff = m.m_tersoff;
	m_m = m.m_m;
	m_I[0] = m.m_I[0];
	m_I[1] = m.m_I[1];
	m_I[2] = m.m_I[2];
	m_invI[0] = m.m_invI[0];
	m_invI[1] = m.m_invI[1];
	m_invI[2] = m.m_invI[2];

	m_numsites = m.m_numsites;
	m_numorientedsites = m.m_numorientedsites;
	assert(m_numsites);
	m_sites_d = new double[m_numsites*3];
	assert(m_sites_d);
	//for(unsigned int i=0;i<m_numsites*3;++i) m_sites_d[i]=m.m_sites_d[i]; // not necessary -> cache only
	m_ljcenters_d = &(m_sites_d[0]);
	m_charges_d = &(m_ljcenters_d[numLJcenters()*3]);
	m_dipoles_d = &(m_charges_d[numCharges()*3]);
	m_quadrupoles_d = &(m_dipoles_d[numDipoles()*3]);
	m_tersoff_d = &(m_quadrupoles_d[numQuadrupoles()*3]);

	m_osites_e = new double[m_numorientedsites*3];
	assert(m_osites_e);
	//for(unsigned int i=0;i<m_numorientedsites*3;++i) m_osites_e[i]=m.m_osites_e[i]; // not necessary -> cache only
	m_dipoles_e = &(m_osites_e[0]);
	m_quadrupoles_e = &(m_dipoles_e[numDipoles()*3]);

	m_sites_F = new double[m_numsites*3];

	assert(m_sites_F);
	//for(unsigned int i=0;i<m_numsites*3;++i) m_sites_F[i]=m.m_sites_F[i]; // not necessary -> cache only
	m_ljcenters_F = &(m_sites_F[0]);
	m_charges_F = &(m_ljcenters_F[numLJcenters()*3]);
	m_dipoles_F = &(m_charges_F[numCharges()*3]);
	m_quadrupoles_F = &(m_dipoles_F[numDipoles()*3]);
	m_tersoff_F = &(m_quadrupoles_F[numQuadrupoles()*3]);
	//m_nextincell=m.m_nextincell;  // not necessary -> temporary only
	m_curTN = 0;
	fixedx = m.fixedx;
	fixedy = m.fixedy;
}

Molecule::Molecule(istream& istrm, streamtype type, const vector<Component>* components) {
	if (type == RESTART) {
		istrm.read((char *) &m_id, sizeof(m_id));
		istrm.read((char *) &m_componentid, sizeof(m_componentid));
		istrm.read((char *) &m_r[0], sizeof(m_r[0]));
		istrm.read((char *) &m_r[1], sizeof(m_r[1]));
		istrm.read((char *) &m_r[2], sizeof(m_r[2]));
		istrm.read((char *) &m_v[0], sizeof(m_v[0]));
		istrm.read((char *) &m_v[1], sizeof(m_v[1]));
		istrm.read((char *) &m_v[2], sizeof(m_v[2]));
		double qw, qx, qy, qz;
		istrm.read((char *) &qw, sizeof(qw));
		istrm.read((char *) &qx, sizeof(qx));
		istrm.read((char *) &qy, sizeof(qy));
		istrm.read((char *) &qz, sizeof(qz));
		m_q = Quaternion(qw, qx, qy, qz);
		istrm.read((char *) &m_D[0], sizeof(m_D[0]));
		istrm.read((char *) &m_D[1], sizeof(m_D[1]));
		istrm.read((char *) &m_D[2], sizeof(m_D[2]));
		istrm.read((char *) &m_F[0], sizeof(m_F[0]));
		istrm.read((char *) &m_F[1], sizeof(m_F[1]));
		istrm.read((char *) &m_F[2], sizeof(m_F[2]));
		istrm.read((char *) &m_M[0], sizeof(m_M[0]));
		istrm.read((char *) &m_M[1], sizeof(m_M[1]));
		istrm.read((char *) &m_M[2], sizeof(m_M[2]));
	}
	m_sites_d = m_sites_F = NULL;
	if (components)
		setupCache(components);
}

double Molecule::dist2(const Molecule& a, double L[3], double dr[]) const {
	double d2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		double L_halve = .5 * L[d];
		dr[d] = a.m_r[d] - m_r[d];
		if (dr[d] > L_halve) {
			dr[d] -= L[d];
		} else if (dr[d] < -L_halve) {
			dr[d] += L[d];
		}
		//dr[d]-=L[d]*floor((dr[d]+L_halve)/L[d]);
		d2 += dr[d] * dr[d];
	}
	return d2;
}

void Molecule::upd_preF(double dt, double vcorr, double Dcorr) {
	assert(m_m);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / m_m;
	for (unsigned short d = 0; d < 3; ++d) {
		m_v[d] = vcorr * m_v[d] + dtInv2m * m_F[d];
		m_oldr[d] = m_r[d];
		m_r[d] += dt * m_v[d];
	}

	double w[3];
	m_q.rotate(m_D, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= m_invI[d];
	Quaternion qhalfstep;
	m_q.differentiate(w, qhalfstep);
	qhalfstep.scale(dt_halve);
	qhalfstep.add(m_q);
	double qcorr = 1. / sqrt(qhalfstep.magnitude2());
	qhalfstep.scale(qcorr);
	for (unsigned short d = 0; d < 3; ++d)
		m_D[d] = Dcorr * m_D[d] + dt_halve * m_M[d];
	qhalfstep.rotate(m_D, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= m_invI[d];
	Quaternion qincr;
	qhalfstep.differentiate(w, qincr);
	qincr.scale(dt);
	m_q.add(qincr);
	qcorr = 1. / sqrt(m_q.magnitude2());
	m_q.scale(qcorr);

}

void Molecule::upd_cache() {
	// update Cache (rotate sites and save relative positions)
	unsigned int i;
	unsigned int ns = numLJcenters();
	for (i = 0; i < ns; ++i)
		m_q.rotateinv((*m_ljcenters)[i].r(), &(m_ljcenters_d[i*3]));
	ns = numCharges();
	for (i = 0; i < ns; ++i)
		m_q.rotateinv((*m_charges)[i].r(), &(m_charges_d[i*3]));
	ns = numDipoles();

	for (i = 0; i < ns; ++i) {
		const Dipole& di = (*m_dipoles)[i];
		m_q.rotateinv(di.r(), &(m_dipoles_d[i*3]));
		m_q.rotateinv(di.e(), &(m_dipoles_e[i*3]));
	}
	ns = numQuadrupoles();

	for (i = 0; i < ns; ++i) {
		const Quadrupole& qi = (*m_quadrupoles)[i];
		m_q.rotateinv(qi.r(), &(m_quadrupoles_d[i*3]));
		m_q.rotateinv(qi.e(), &(m_quadrupoles_e[i*3]));
	}
	ns = this->numTersoff();
	for (i = 0; i < ns; i++)
		m_q.rotateinv((*m_tersoff)[i].r(), &(m_tersoff_d[i*3]));

	clearFM();
	//m_Upot=0.;
}


void Molecule::upd_postF(double dt_halve, double& summv2, double& sumIw2) {
	//m_Upot*=.5;


	calcFM();

	//if(m_id==1) cout << "Kraft: " << m_F[0] << " / " << m_F[1] << " / " << m_F[2] << endl;

	double dtInv2m = dt_halve / m_m;
	double v2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		m_v[d] += dtInv2m * m_F[d];
		v2 += m_v[d] * m_v[d];
		m_D[d] += dt_halve * m_M[d];
	}
#ifndef NDEBUG
	if (!(v2 > 0.0)) {
		/*
		 * catches NaN forces and coordinates
		 */
		global_log->error() << "Molecule " << m_id << " has v^2 = " << v2 << ".\n"
		                    << "r = (" << m_r[0] << " / " << m_r[1] << " / " << m_r[2] << "), "
		                    << "F = (" << m_F[0] << " / " << m_F[1] << " / " << m_F[2] << "), "
		                    << "v = (" << m_v[0] << " / " << m_v[1] << " / " << m_v[2] << ")." << endl;
		exit(1);
	}
#endif
	summv2 += m_m * v2;
	double w[3];
	m_q.rotate(m_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= m_invI[d];
		Iw2 += m_I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

double Molecule::Urot() {
	double w[3];
	m_q.rotate(m_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= m_invI[d];
		Iw2 += m_I[d] * w[d] * w[d];
	}
	return 0.5 * Iw2;
}

void Molecule::calculate_mv2_Iw2(double& summv2, double& sumIw2) {
	summv2 += m_m * v2();
	double w[3];
	m_q.rotate(m_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= m_invI[d];
		Iw2 += m_I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void Molecule::calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) {
	double vcx = m_v[0] - offx;
	double vcy = m_v[1] - offy;
	double vcz = m_v[2] - offz;
	summv2 += m_m * (vcx*vcx + vcy*vcy + vcz*vcz);

	double w[3];
	m_q.rotate(m_D, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= m_invI[d];
		Iw2 += m_I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void Molecule::scale_v(double s, double offx, double offy, double offz) {
	this->vsub(offx, offy, offz);
	this->scale_v(s);
	this->vadd(offx, offy, offz);
}

void Molecule::write(ostream& ostrm) const {
	ostrm << m_id << "\t" << (m_componentid+1) << "\t"
	      << m_r[0] << " " << m_r[1] << " " << m_r[2] << "\t"
	      << m_v[0] << " " << m_v[1] << " " << m_v[2] << "\t"
	      << m_q.qw() << " " << m_q.qx() << " " << m_q.qy() << " " << m_q.qz() << "\t"
	      << m_D[0] << " " << m_D[1] << " " << m_D[2] << "\t"
	      << endl;
}

void Molecule::save_restart(ostream& ostrm) const {
	ostrm.write((const char *) &m_id, sizeof(m_id));
	ostrm.write((const char *) &m_componentid, sizeof(m_componentid));
	ostrm.write((const char *) &m_r[0], sizeof(m_r[0]));
	ostrm.write((const char *) &m_r[1], sizeof(m_r[1]));
	ostrm.write((const char *) &m_r[2], sizeof(m_r[2]));
	ostrm.write((const char *) &m_v[0], sizeof(m_v[0]));
	ostrm.write((const char *) &m_v[1], sizeof(m_v[1]));
	ostrm.write((const char *) &m_v[2], sizeof(m_v[2]));
	double q = m_q.qw();
	ostrm.write((const char *) &q, sizeof(q));
	q = m_q.qx();
	ostrm.write((const char *) &q, sizeof(q));
	q = m_q.qy();
	ostrm.write((const char *) &q, sizeof(q));
	q = m_q.qz();
	ostrm.write((const char *) &q, sizeof(q));
	ostrm.write((const char *) &m_D[0], sizeof(m_D[0]));
	ostrm.write((const char *) &m_D[1], sizeof(m_D[1]));
	ostrm.write((const char *) &m_D[2], sizeof(m_D[2]));
	ostrm.write((const char *) &m_F[0], sizeof(m_F[0]));
	ostrm.write((const char *) &m_F[1], sizeof(m_F[1]));
	ostrm.write((const char *) &m_F[2], sizeof(m_F[2]));
	ostrm.write((const char *) &m_M[0], sizeof(m_M[0]));
	ostrm.write((const char *) &m_M[1], sizeof(m_M[1]));
	ostrm.write((const char *) &m_M[2], sizeof(m_M[2]));
}

void Molecule::addTersoffNeighbour(Molecule* m, bool pairType) {
	// this->m_Tersoff_neighbours.insert(pair<Molecule*, bool>(m, (pairType > 0)));
	for (int j = 0; j < m_curTN; j++) {
		if (m->m_id == m_Tersoff_neighbours_first[j]->id()) {
			this->m_Tersoff_neighbours_first[j] = m;
			this->m_Tersoff_neighbours_second[j] = pairType;
			return;
		}
	}

	this->m_Tersoff_neighbours_first[m_curTN] = m;
	this->m_Tersoff_neighbours_second[m_curTN] = pairType;
	this->m_curTN++;
	if (m_curTN > MAXTN) {
		global_log->error() << "Tersoff neighbour list overflow: Molecule " << m->m_id << " has more than " << MAXTN << " Tersoff neighbours." << endl;
		exit(1);
	}
}

double Molecule::tersoffParameters(double params[15]) //returns delta_r
{
	const Tersoff* t = &((*(this->m_tersoff))[0]);
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
	params[10] = 3.141592653589 / (t->S() - t->R());
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
	assert(m_componentid >= 0);
	m_numsites = m_numorientedsites = 0;
	m_ljcenters = &(*components)[m_componentid].ljcenters();
	m_numsites += m_ljcenters->size();
	m_charges = &(*components)[m_componentid].charges();
	m_numsites += m_charges->size();
	m_dipoles = &(*components)[m_componentid].dipoles();
	m_numsites += m_dipoles->size();
	m_numorientedsites += m_dipoles->size();
	m_quadrupoles = &(*components)[m_componentid].quadrupoles();
	m_numsites += m_quadrupoles->size();
	m_numorientedsites += m_quadrupoles->size();
	m_tersoff = &(*components)[m_componentid].tersoff();
#ifndef NDEBUG
	if (!m_tersoff) {
		global_log->error() << "Tersoff vector null pointer detected for Molecule " << m_id << endl;
		exit(1);
	}
#endif
	m_numsites += m_tersoff->size();

	m_m = (*components)[m_componentid].m();
	m_I[0] = (*components)[m_componentid].I11();
	m_I[1] = (*components)[m_componentid].I22();
	m_I[2] = (*components)[m_componentid].I33();
	for (unsigned short d = 0; d < 3; ++d) {
		if (m_I[d] != 0.)
			m_invI[d] = 1. / m_I[d];
		else
			m_invI[d] = 0.;
	}

	assert(m_numsites);

	m_sites_d = new double[m_numsites*3];
	assert(m_sites_d);
	m_ljcenters_d = &(m_sites_d[0]);
	m_charges_d = &(m_ljcenters_d[numLJcenters()*3]);
	m_dipoles_d = &(m_charges_d[numCharges()*3]);
	m_quadrupoles_d = &(m_dipoles_d[numDipoles()*3]);
	m_tersoff_d = &(m_quadrupoles_d[numQuadrupoles()*3]);

	m_osites_e = new double[m_numorientedsites*3];
	assert(m_osites_e);
	m_dipoles_e = &(m_osites_e[0]);
	m_quadrupoles_e = &(m_dipoles_e[numDipoles()*3]);

	m_sites_F = new double[m_numsites*3];

	assert(m_sites_F);
	m_ljcenters_F = &(m_sites_F[0]);
	m_charges_F = &(m_ljcenters_F[numLJcenters()*3]);
	m_dipoles_F = &(m_charges_F[numCharges()*3]);
	m_quadrupoles_F = &(m_dipoles_F[numDipoles()*3]);
	m_tersoff_F = &(m_quadrupoles_F[numQuadrupoles()*3]);

	this->clearFM();
}

inline void Molecule::clearFM() {
	for (unsigned int i = 0; i < m_numsites * 3; ++i)
		m_sites_F[i] = 0.;
	m_F[0] = m_F[1] = m_F[2] = 0.;
	m_M[0] = m_M[1] = m_M[2] = 0.;
}

inline void Molecule::calcFM() {
	unsigned int ns = numSites();
	for (unsigned int si = 0; si < ns; ++si) {
		const double* Fsite = site_F(si);
		const double* dsite = site_d(si);
#ifndef NDEBUG
		/*
		 * catches NaN assignments
		 */
		for (int d = 0; d < 3; d++) {
			if (!((dsite[d] >= 0) || (dsite[d] < 0))) {
				global_log->error() << "Severe dsite[" << d << "] error for site " << si << " of m" << m_id << endl;
				assert(false);
			}
			if (!((Fsite[d] >= 0) || (Fsite[d] < 0))) {
				global_log->error() << "Severe Fsite[" << d << "] error for site " << si << " of m" << m_id << endl;
				assert(false);
			}
		}
#endif
		Fadd(Fsite);
		m_M[0] += dsite[1] * Fsite[2] - dsite[2] * Fsite[1];
		m_M[1] += dsite[2] * Fsite[0] - dsite[0] * Fsite[2];
		m_M[2] += dsite[0] * Fsite[1] - dsite[1] * Fsite[0];
	}
}

void Molecule::setDomain(Domain* domain) {
	Molecule::_domain = domain;
}

/*
 * catches NaN values and missing data
 */
void Molecule::check(unsigned long id) {
	assert(m_id == id);
	assert(m_componentid >= 0);
	assert(m_m > 0.0);
	assert(m_numsites > 0);
	assert(m_numorientedsites >= 0);
	for (int d = 0; d < 3; d++) {
		assert((m_r[d] >= 0.0) || (m_r[d] < 0.0));
		assert((m_v[d] >= 0.0) || (m_v[d] < 0.0));
		assert((m_D[d] >= 0.0) || (m_D[d] < 0.0));
		assert((m_F[d] >= 0.0) || (m_F[d] < 0.0));
		assert((m_M[d] >= 0.0) || (m_M[d] < 0.0));
		assert((m_I[d] >= 0.0) || (m_I[d] < 0.0));
		assert((m_invI[d] >= 0.0) || (m_invI[d] < 0.0));
	}
}
