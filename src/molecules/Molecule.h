/*************************************************************************
 * Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or (at *
 * your option) any later version.                                       *
 *                                                                       *
 * This program is distributed in the hope that it will be useful, but   *
 * WITHOUT ANY WARRANTY; without even the implied warranty of            * 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      *
 * General Public License for more details.                              *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the Free Software           *
 * Foundation, 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.   *
 *************************************************************************/

#ifndef MOLECULE_H_
#define MOLECULE_H_

/*
 * maximal size of the Tersoff neighbour list
 */
#define MAXTN 10

#include "molecules/Quaternion.h"
#include "molecules/Comp2Param.h"
#include "molecules/Site.h"
#include "molecules/Component.h"
#include "integrators/Integrator.h"
class Domain;

#include <vector>
#include <iostream>
#include <string>

#include <cassert>

//! @brief Molecule modeled as LJ sphere with point polarities + Tersoff potential
//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
class Molecule {

public:
	enum streamtype {
		RESTART
	};

	Molecule(unsigned long id = 0, int componentid = 0,
	         double rx = 0., double ry = 0., double rz = 0.,
	         double vx = 0., double vy = 0., double vz = 0.,
	         double q0 = 0., double q1 = 0., double q2 = 0., double q3 = 0.,
	         double Dx = 0., double Dy = 0., double Dz = 0.,
	         const std::vector<Component>* components = NULL
	);
	Molecule(const Molecule& m);
	Molecule(std::istream& istrm, streamtype type, const std::vector<Component>* components=NULL);

	~Molecule() {
		assert(_sites_d); delete[] _sites_d;
		assert(_osites_e); delete[] _osites_e;
		assert(_sites_F); delete[] _sites_F;
	}

	/** get the ID */
	unsigned long id() const { return _id; }
	void setid(unsigned long id) { this->_id = id; }
	/** get the Component */
	int componentid() const { return _componentid; }
	/** get the position */
	double r(unsigned short d) const { return _r[d]; }

	double oldr (unsigned short d) const { return _oldr[d]; }

	/** get the velocity */
	double v(unsigned short d) const { return _v[d]; }
	/** get the Orientation */
	const Quaternion& q() const { return _q; }

	inline void move(int d, double dr) { _r[d] += dr; }

	/** get the rotatational speed */
	double D(unsigned short d) const { return _D[d]; }

	/** get F */
	double F(unsigned short d) const {return _F[d]; }
	/** get M */
	double M(unsigned short d) const {return _M[d]; }

	//double Upot() const { return _Upot; }
	double Utrans() const { return .5*_m*(_v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2]); }
	double Urot();

	/** get number of sites */
	unsigned int numSites() const { return _numsites; }
	unsigned int numLJcenters() const { return _ljcenters->size(); }
	unsigned int numCharges() const { return _charges->size(); }
	unsigned int numDipoles() const { return _dipoles->size(); }
	unsigned int numQuadrupoles() const { return _quadrupoles->size(); }
	unsigned int numTersoff() const { assert(_tersoff); return _tersoff->size(); }

	const double* site_d(unsigned int i) const { return &(_sites_d[3*i]); }
	const double* osite_e(unsigned int i) const { return &(_osites_e[3*i]); }
	const double* site_F(unsigned int i) const { return &(_sites_F[3*i]); }
	const double* ljcenter_d(unsigned int i) const { return &(_ljcenters_d[3*i]); }
	const double* ljcenter_F(unsigned int i) const { return &(_ljcenters_F[3*i]); }
	const double* charge_d(unsigned int i) const { return &(_charges_d[3*i]); }
	const double* charge_F(unsigned int i) const { return &(_charges_F[3*i]); }
	const double* dipole_d(unsigned int i) const { return &(_dipoles_d[3*i]); }
	const double* dipole_e(unsigned int i) const { return &(_dipoles_e[3*i]); }
	const double* dipole_F(unsigned int i) const { return &(_dipoles_F[3*i]); }
	const double* quadrupole_d(unsigned int i) const { return &(_quadrupoles_d[3*i]); }
	const double* quadrupole_e(unsigned int i) const { return &(_quadrupoles_e[3*i]); }
	const double* quadrupole_F(unsigned int i) const { return &(_quadrupoles_F[3*i]); }
	const double* tersoff_d(unsigned int i) const { return &(_tersoff_d[3*i]); }
	const double* tersoff_F(unsigned int i) const { return &(_tersoff_F[3*i]); }

	/** get object memory size */
	static unsigned long memsize() { return sizeof(Molecule); }

	/** set the position */
	void setr(unsigned short d, double r) { _r[d]=r; }

	void setOldRFromR () {
		for (int i = 0; i < 3; i++)
			_oldr[i] = _r[i];
	}

	/** calculate the difference vector and return the square (euclidean) distance */
	double dist2(const Molecule& a, double dr[]) const {
		double d2=0.;
		for (unsigned short d=0; d<3; ++d) {
			dr[d] = a._r[d] - _r[d];
			d2 += dr[d] * dr[d];
		}
		return d2;
	}
	double dist2(const Molecule& a, double L[3], double dr[]) const;
	/** calculate and return the square velocity */
	double v2() const {return _v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2]; }

	/** set molecule force and momentum 
	 * @param F force vector (x,y,z)
	 * @param M momentum vector (x,y,z)
	 */
	void setF(double F[3]) { for(int d = 0; d < 3; d++ ) { _F[d] = F[d]; } }
	void setM(double M[3]) { for(int d = 0; d < 3; d++ ) { _M[d] = M[d]; } }

	void scale_v(double s) { for(unsigned short d=0;d<3;++d) _v[d]*=s; }
	void scale_v(double s, double offx, double offy, double offz);
	void scale_F(double s) { for(unsigned short d=0;d<3;++d) _F[d]*=s; }
	void scale_D(double s) { for(unsigned short d=0;d<3;++d) _D[d]*=s; }
	void scale_M(double s) { for(unsigned short d=0;d<3;++d) _M[d]*=s; }

	void Fadd(const double a[]) { for(unsigned short d=0;d<3;++d) _F[d]+=a[d]; }
	void Fsub(const double a[]) { for(unsigned short d=0;d<3;++d) _F[d]-=a[d]; }

	void Madd(const double a[]) { for(unsigned short d=0;d<3;++d) _M[d]+=a[d]; }
	void Msub(const double a[]) { for(unsigned short d=0;d<3;++d) _M[d]-=a[d]; }

	void vadd(const double ax, const double ay, const double az) {
		_v[0] += ax; _v[1] += ay; _v[2] += az;
	}
	void vsub(const double ax, const double ay, const double az) {
		_v[0] -= ax; _v[1] -= ay; _v[2] -= az;
	}
	void setXY() { fixedx = _r[0]; fixedy = _r[1]; }
	void resetXY()
	{
		_v[0] = 0.0;
		_v[1] = 0.0;
		_F[1] = 0.0;
		_F[0] = 0.0;
		_r[0] = fixedx;
		_r[1] = fixedy;
	}

	void Fsiteadd(unsigned int i, double a[])
	{ double* Fsite=&(_sites_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fsitesub(unsigned int i, double a[])
	{ double* Fsite=&(_sites_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Fljcenteradd(unsigned int i, double a[])
	{ double* Fsite=&(_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fljcentersub(unsigned int i, double a[])
	{ double* Fsite=&(_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Fchargeadd(unsigned int i, double a[])
	{ double* Fsite=&(_charges_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fchargesub(unsigned int i, double a[])
	{ double* Fsite=&(_charges_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Fdipoleadd(unsigned int i, double a[])
	{ double* Fsite=&(_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fdipolesub(unsigned int i, double a[])
	{ double* Fsite=&(_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Fquadrupoleadd(unsigned int i, double a[])
	{ double* Fsite=&(_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Fquadrupolesub(unsigned int i, double a[])
	{ double* Fsite=&(_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
	void Ftersoffadd(unsigned int i, double a[])
	{ double* Fsite=&(_tersoff_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
	void Ftersoffsub(unsigned int i, double a[])
	{ double* Fsite=&(_tersoff_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }

	void upd_preF(double dt, double vcorr=1., double Dcorr=1.);
	void upd_cache();
	void upd_postF(double dt_halve, double& summv2, double& sumIw2);

	//! calculate summv2 and sumIw2
	//! @todo what is sumIw2?
	                         //! @todo comment
	void calculate_mv2_Iw2(double& summv2, double& sumIw2);
	void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz);

	/** write information to stream */
	void write(std::ostream& ostrm) const;
	/** write binary information to stream */
	void save_restart(std::ostream& ostrm) const;

	static void setDomain(Domain* domain);

	inline unsigned getCurTN() { return this->_curTN; }
	inline Molecule* getTersoffNeighbour(unsigned i) { return this->_Tersoff_neighbours_first[i]; }
	inline bool getPairCode(unsigned i) { return this->_Tersoff_neighbours_second[i]; }
	inline void clearTersoffNeighbourList() { this->_curTN = 0; }
	void addTersoffNeighbour(Molecule* m, bool pairType);
	double tersoffParameters(double params[15]); //returns delta_r

	// clear forces and moments
	void clearFM();
	void check(unsigned long id);

private:

	static Domain* _domain;
	unsigned long _id; // IDentification number of that molecule
	int _componentid;  // IDentification number of its component type
	double _r[3];  // position coordinates
	double _oldr[3]; // position coordinates last step
	double _v[3];  // velocity
	Quaternion _q; // orientation
	double _D[3];  // angular momentum

	double _F[3];  // forces
	double _M[3];  // moments

	const std::vector<LJcenter>* _ljcenters;
	const std::vector<Charge>* _charges;
	const std::vector<Dipole>* _dipoles;
	const std::vector<Quadrupole>* _quadrupoles;
	const std::vector<Tersoff>* _tersoff;

	double _m; // total mass
	double _I[3],_invI[3];  // moment of inertia for principal axes and it's inverse
	std::size_t _numsites; // number of sites
	std::size_t _numorientedsites; // number of oriented sites (subset of sites)
	// global site coordinates relative to site origin
	// row order: dx1,dy1,dz1,dx2,dy2,dz2,...
	double *_sites_d;
	double *_ljcenters_d, *_charges_d, *_dipoles_d,
	       *_quadrupoles_d, *_tersoff_d;
	// site orientation
	double *_osites_e;
	double *_dipoles_e, *_quadrupoles_e;
	// site Forces
	// row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
	double* _sites_F;
	double *_ljcenters_F, *_charges_F, *_dipoles_F,
	       *_quadrupoles_F, *_tersoff_F;

	Molecule* _Tersoff_neighbours_first[MAXTN];
	bool _Tersoff_neighbours_second[MAXTN];
	int _curTN;
	double fixedx, fixedy;

	// setup cache values/properties
	void setupCache(const std::vector<Component>* components);
	// calculate forces and moments for already given site forces
	void calcFM();
};

#endif /*MOLECULE_H_*/
