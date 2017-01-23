#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <vector>
#include <iostream>
#include <cassert>
#include <array>
#include <string>

#include "molecules/Component.h"
#include "molecules/Comp2Param.h"
#include "molecules/Quaternion.h"
#include "molecules/Site.h"


class Domain;
class CellDataSoA;

//! @brief Molecule modeled as LJ sphere with point polarities
class Molecule {

public:
	// TODO Correct this constructor: the components vector is optional,
	// but if it is left away, all pointer data is not initialized (which is not
	// neccessarily bad), but then assertions fail (e.g. in the destructor) and we can't
	// use it's instances.
	Molecule(unsigned long id = 0, Component *component = nullptr,
	         double rx = 0., double ry = 0., double rz = 0.,
	         double vx = 0., double vy = 0., double vz = 0.,
	         double q0 = 0., double q1 = 0., double q2 = 0., double q3 = 0.,
	         double Dx = 0., double Dy = 0., double Dz = 0.
	);
	Molecule(const Molecule& m);

	Molecule& operator=(const Molecule& m);

	~Molecule() {
		// don't delete SoA
		_soa = nullptr;
	}

	/** get molecule ID */
	unsigned long id() const { return _id; }
	/** set molecule ID */
	void setid(unsigned long id) { _id = id; }
	/** get the molecule's component ID */
	unsigned int componentid() const { return _component->ID(); }
	/** set the molecule's component */
	void setComponent(Component *component) { _component = component; this->updateMassInertia();}
	/** return pointer to component to which the molecule belongs */
	Component* component() const { return _component; }
	/** get component lookUpID */
	unsigned getComponentLookUpID() const { return _component->getLookUpId();}
	/** get position coordinate */
	double r(unsigned short d) const { return _r[d]; }
	/** set position coordinate */
	void setr(unsigned short d, double r) { _r[d] = r; }
	/** get velocity coordinate */
	double v(unsigned short d) const { return _v[d]; }
	/** set velocity */
	void setv(unsigned short d, double v) { _v[d] = v; }
	/** get molecule's mass */
	double mass() const { return _m; }

	/** get coordinate of current force onto molecule */
	double F(unsigned short d) const {return _F[d]; }
	/** get molecule's orientation */
	const Quaternion& q() const { return _q; }


	/** set molecule's orientation */
	void setq(Quaternion q){ _q = q; }

	/** get coordinate of the rotatational speed */
	double D(unsigned short d) const { return _L[d]; }
	/** get coordinate of the current angular momentum  onto molecule */ 
	double M(unsigned short d) const { return _M[d]; }
	/** get the virial **/
	double Vi(unsigned short d) const { return _Vi[d];}

	void setD(unsigned short d, double D) { this->_L[d] = D; }

	inline void move(int d, double dr) { _r[d] += dr; }

	// by Stefan Becker <stefan.becker@mv.uni-kl.de> 
	// method returns the total mass of a particle
	double gMass(){return _m;}
	//by Stefan Becker
		/** get the moment of inertia of a particle */
	double getI(unsigned short d) const { return _I[d]; }
	/** update mass and moment of inertia by component definition */
	void updateMassInertia() {
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
	}

	/** calculate and return the square velocity */
	double v2() const {return _v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2]; }
	
	/** return the translational energy of the molecule */
	double U_trans() const { return 0.5 * _m * v2(); }
	/** return the rotational energy of the molecule */
	double U_rot();
	/** return total kinetic energy of the molecule */
	double U_kin() { return U_trans() + U_rot(); }
	
	void setupSoACache(CellDataSoA * const s, unsigned iLJ, unsigned iC, unsigned iD, unsigned iQ);

	void setSoA(CellDataSoA * const s) {_soa = s;}
	void setStartIndexSoA_LJ(unsigned i) {_soa_index_lj = i;}
	void setStartIndexSoA_C(unsigned i) {_soa_index_c = i;}
	void setStartIndexSoA_D(unsigned i) {_soa_index_d = i;}
	void setStartIndexSoA_Q(unsigned i) {_soa_index_q = i;}

	/* TODO: Maybe we should better do this using the component directly? 
	 * In the GNU STL vector.size() causes two memory accesses and one subtraction!
	 */
	/** get number of sites */
	unsigned int numSites() const { return _component->numSites(); }
	unsigned int numOrientedSites() const { return _component->numOrientedSites();  }
	unsigned int numLJcenters() const { return _component->numLJcenters(); }
	unsigned int numCharges() const { return _component->numCharges(); }
	unsigned int numDipoles() const { return _component->numDipoles(); }
	unsigned int numQuadrupoles() const { return _component->numQuadrupoles(); }

	std::array<double, 3> site_d(unsigned int i) const {
		const unsigned n1 = numLJcenters(), n2 = numCharges()+ n1, n3 = numDipoles() + n2;
#ifndef NDEBUG
		const unsigned n4 = numQuadrupoles() + n3;
#endif
		if(i < n1) {
			return ljcenter_d(i);
		} else if (i < n2) {
			return charge_d(i - n1);
		} else if (i < n3) {
			return dipole_d(i - n2);
		} else { assert(i < n4);
			return quadrupole_d(i - n3);
		}
	}
	std::array<double, 3> ljcenter_d(unsigned int i) const {
		std::array<double, 3> ret;
		computeLJcenter_d(i,ret.data());
		return ret;
	}
	std::array<double, 3> charge_d(unsigned int i) const {
		std::array<double, 3> ret;
		computeCharge_d(i,ret.data());
		return ret;
	}
	std::array<double, 3> dipole_d(unsigned int i) const {
		std::array<double, 3> ret;
		computeDipole_d(i, ret.data());
		return ret;
	}
	std::array<double, 3> quadrupole_d(unsigned int i) const {
		std::array<double, 3> ret;
		computeQuadrupole_d(i, ret.data());
		return ret;
	}

	std::array<double, 3> site_d_abs(unsigned int i) const {
		const unsigned n1 = numLJcenters(), n2 = numCharges()+ n1, n3 = numDipoles() + n2;
#ifndef NDEBUG
		const unsigned n4 = numQuadrupoles() + n3;
#endif
		if(i < n1) {
			return ljcenter_d_abs(i);
		} else if (i < n2) {
			return charge_d_abs(i - n1);
		} else if (i < n3) {
			return dipole_d_abs(i - n2);
		} else { assert(i < n4);
			return quadrupole_d_abs(i - n3);
		}
	}
	std::array<double, 3> ljcenter_d_abs(unsigned int i) const;
	std::array<double, 3> charge_d_abs(unsigned int i) const;
	std::array<double, 3> dipole_d_abs(unsigned int i) const;
	std::array<double, 3> quadrupole_d_abs(unsigned int i) const;

	std::array<double, 3> dipole_e(unsigned int i) const;
	std::array<double, 3> quadrupole_e(unsigned int i) const;

	std::array<double, 3> site_F(unsigned int i) const {
		const unsigned n1 = numLJcenters(), n2 = numCharges()+ n1, n3 = numDipoles() + n2;
#ifndef NDEBUG
		const unsigned n4 = numQuadrupoles() + n3;
#endif
		if(i < n1) {
			return ljcenter_F(i);
		} else if (i < n2) {
			return charge_F(i - n1);
		} else if (i < n3) {
			return dipole_F(i - n2);
		} else { assert(i < n4);
			return quadrupole_F(i - n3);
		}
	}
	std::array<double, 3> ljcenter_F(unsigned int i) const;
	std::array<double, 3> charge_F(unsigned int i) const;
	std::array<double, 3> dipole_F(unsigned int i) const;
	std::array<double, 3> quadrupole_F(unsigned int i) const;

	void normalizeQuaternion() {
		_q.normalize();
	}
	void computeLJcenter_d(unsigned int i, double result[3]) const {
		assert(_q.isNormalized());
		_q.rotate(_component->ljcenter(i).r(), result);
	}
	void computeCharge_d(unsigned int i, double result[3]) const {
		assert(_q.isNormalized());
		_q.rotate(_component->charge(i).r(), result);
	}
	void computeDipole_d(unsigned int i, double result[3]) const {
		assert(_q.isNormalized());
		_q.rotate(_component->dipole(i).r(), result);
	}
	void computeQuadrupole_d(unsigned int i, double result[3]) const {
		assert(_q.isNormalized());
		_q.rotate(_component->quadrupole(i).r(), result);
	}
	void computeDipole_e(unsigned int i, double result[3]) const {
		assert(_q.isNormalized());
		_q.rotate(_component->dipole(i).e(), result);
	}
	void computeQuadrupole_e(unsigned int i, double result[3]) const {
		assert(_q.isNormalized());
		_q.rotate(_component->quadrupole(i).e(), result);
	}


	/**
	 * get the total object memory size, together with all its members
	 * \Note You can retrieve the size of the molecule class itself simply
	 *       with the sizeof()-operator.
	 */
	unsigned long totalMemsize() const;


	/* TODO: Is this realy necessary? Better use a function like dist2(m1.r(), m2.r()). */
	/** Calculate the difference vector and return the square (euclidean) distance.
	 *
	 *  \param molecule2 molecule to which the distance shall be calculated
	 */
	double dist2(const Molecule& molecule2, double dr[3]) const {
		double d2 = 0.;
		for (unsigned short d = 0; d < 3; d++) {
			dr[d] = molecule2._r[d] - _r[d];
			d2 += dr[d] * dr[d];
		}
		return d2;
	}
	

	/** set force acting on molecule
	 * @param[out] F force vector (x,y,z)
	 */
	void setF(double F[3]) { for(int d = 0; d < 3; d++ ) { _F[d] = F[d]; } }

	/** set momentum acting on molecule 
	 * @param M force vector (x,y,z)
	 */
	void setM(double M[3]) { for(int d = 0; d < 3; d++ ) { _M[d] = M[d]; } }
	void setVi(double Vi[3]) { for(int d = 0; d < 3; d++) { _Vi[d] = Vi[d]; } }

	void scale_v(double s) { for(unsigned short d=0;d<3;++d) _v[d]*=s; }
	void scale_v(double s, double offx, double offy, double offz);
	void scale_F(double s) { for(unsigned short d=0;d<3;++d) _F[d]*=s; }
	void scale_D(double s) { for(unsigned short d=0;d<3;++d) _L[d]*=s; }
	void scale_M(double s) { for(unsigned short d=0;d<3;++d) _M[d]*=s; }

	void Fadd(const double a[]) { for(unsigned short d=0;d<3;++d) _F[d]+=a[d]; }

	void Madd(const double a[]) { for(unsigned short d=0;d<3;++d) _M[d]+=a[d]; }
	
	void Viadd(const double a[]) { for(unsigned short d=0;d<3;++d) _Vi[d]+=a[d]; }

	void vadd(const double ax, const double ay, const double az) {
		_v[0] += ax; _v[1] += ay; _v[2] += az;
	}
	void vsub(const double ax, const double ay, const double az) {
		_v[0] -= ax; _v[1] -= ay; _v[2] -= az;
	}

	void Fljcenteradd(unsigned int i, double a[]);
	void Fljcentersub(unsigned int i, double a[]);
	void Fchargeadd(unsigned int i, double a[]);
	void Fchargesub(unsigned int i, double a[]);
	void Fdipoleadd(unsigned int i, double a[]);
	void Fdipolesub(unsigned int i, double a[]);
	void Fquadrupoleadd(unsigned int i, double a[]);
	void Fquadrupolesub(unsigned int i, double a[]);

	/** First step of the leap frog integrator */
	void upd_preF(double dt);
	/** second step of the leap frog integrator */
	void upd_postF(double dt_halve, double& summv2, double& sumIw2);

	/** @brief Calculate twice the translational and rotational kinetic energies
	 * @param[out] summv2   twice the translational kinetic energy \f$ m v^2 \f$
	 * @param[out] sumIw2   twice the rotational kinetic energy \f$ I \omega^2 \f$
	 */
	void calculate_mv2_Iw2(double& summv2, double& sumIw2);
	void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz);

	/**
	 * @return format of the function write(...)
	 */
	static std::string getWriteFormat();

	/** write information to stream */
	void write(std::ostream& ostrm) const;

	/** clear forces and moments */
	void clearFM();
	/** calculate forces and moments for already given site forces */
	void calcFM();
	
	/** perform data consistency check for the molecule (only debug mode) */
	void check(unsigned long id);

	//! @brief find out whether m1 is before m2 (in some global ordering)
	//!
	//! Compares this molecule to m2 based on their coordinates.
	//!
	//! @return true if this molecule is smaller than m2 (according to the order
	//!         induced by the coordinates (z,y,x)
	//!
	//! At the boundary between two processes (if used in parallel mode), the forces
	//! for pairs which cross the boundary are calculated twice (once by each proc who
	//! owns one of the particles). But the contribution to macroscopic value must be
	//! counted only once, which is done by the process who owns the "first" particle.
	//! As order criterion, the spacial position is used int this method. The particles
	//! with lower x-coordinate is first (if equal, then y- or z-coordinate).
	//! For pairs which are completely on one process, the first particle can be
	//! determined from the cell structure. But for pairs on different procs, the
	//! corresponding cell discretisations might be different as well, and therefore
	//! the cell structure must not be used to determine the order.
	bool isLessThan(const Molecule& m2) const;

	/**
	 * \brief test whether molecule is inside a cuboid region
	 * @param l lower left front corner of cube (equality allowed)
	 * @param u upper right back corner of cube (equality not allowed)
	 * @return true if molecule is contained in the box, false otherwise
	 */
	bool inBox(const double l[3], const double u[3]) const
	{bool in = true; for(int d=0; d < 3; ++d) {in &= (_r[d] >= l[d] and _r[d] < u[d]);} return in;}

private:
    Component *_component;  /**< IDentification number of its component type */
	double _r[3];  /**< position coordinates */
	double _F[3];  /**< forces */
	double _v[3];  /**< velocity */
	Quaternion _q; /**< angular orientation */
	double _M[3];  /**< torsional moment */
	double _L[3];  /**< angular momentum */
	double _Vi[3]; /** Virial tensor **/
    unsigned long _id;  /**< IDentification number of that molecule */

	double _m; /**< total mass */
	double _I[3],_invI[3];  // moment of inertia for principal axes and it's inverse

	/* absolute positions are stored in the soa. Work-arounds to get the relative ones*/
	CellDataSoA * _soa;
	unsigned _soa_index_lj;
	unsigned _soa_index_c;
	unsigned _soa_index_d;
	unsigned _soa_index_q;
};


std::ostream& operator<<( std::ostream& os, const Molecule& m );



/** @brief Calculate the distance between two sites of two molecules.
 *
 * @param[in]  drm distance vector between the two molecule centers
 * @param[in]  ds1 distance vector from the center of molecule1 to its site
 * @param[in]  ds2 distance vector from the center of molecule2 to its site
 * @param[out] drs distance vector site-site
 * @param[out] dr2 distance site-site
 *
 */
inline void SiteSiteDistance(const double drm[3], const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
	for (unsigned short d = 0; d < 3; ++d)
		drs[d] = drm[d] + ds1[d] - ds2[d];
	dr2 = drs[0]*drs[0] + drs[1]*drs[1] + drs[2]*drs[2];
}

/** @brief Calculate the distance between two sites of two molecules.
 *
 * @param[in]  ds1 absolute position of site 1
 * @param[in]  ds2 absolute position of site 2
 * @param[out] drs distance vector site-site
 * @param[out] dr2 distance site-site
 *
 */
inline void SiteSiteDistanceAbs(const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
	for (unsigned short d = 0; d < 3; ++d)
		drs[d] = ds1[d] - ds2[d];
	dr2 = drs[0]*drs[0] + drs[1]*drs[1] + drs[2]*drs[2];
}


inline void minusSiteSiteDistance(const double drm[3], const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
	for (unsigned short d = 0; d < 3; ++d)
		drs[d] = ds2[d] - drm[d] - ds1[d];
	dr2 = drs[0]*drs[0] + drs[1]*drs[1] + drs[2]*drs[2];
}

inline void minusSiteSiteDistanceAbs(const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
	for (unsigned short d = 0; d < 3; ++d)
		drs[d] = ds2[d] - ds1[d];
	dr2 = drs[0]*drs[0] + drs[1]*drs[1] + drs[2]*drs[2];
}

#endif /* MOLECULE_H_ */
