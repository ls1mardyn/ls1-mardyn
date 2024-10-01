#ifndef FULLMOLECULE_H_
#define FULLMOLECULE_H_

#include "MoleculeInterface.h"
#include "molecules/Comp2Param.h"
#include "molecules/Site.h"

#include <vector>
#include <iostream>
#include "utils/mardyn_assert.h"
#include <string>


class Domain;
class CellDataSoA;

//! @brief FullMolecule modeled as LJ sphere with point polarities
class FullMolecule : public MoleculeInterface {

public:
	// TODO Correct this constructor: the components vector is optional,
	// but if it is left away, all pointer data is not initialized (which is not
	// neccessarily bad), but then assertions fail (e.g. in the destructor) and we can't
	// use it's instances.
	FullMolecule(unsigned long id = 0, Component *component = nullptr,
	         double rx = 0., double ry = 0., double rz = 0.,
	         double vx = 0., double vy = 0., double vz = 0.,
	         double q0 = 1., double q1 = 1., double q2 = 0., double q3 = 0.,
	         double Dx = 0., double Dy = 0., double Dz = 0.
	);
	FullMolecule(const FullMolecule& m);

	FullMolecule& operator=(const FullMolecule& m);

	~FullMolecule() override {
		// don't delete SoA
		_soa = nullptr;
	}

	/** get molecule ID */
	unsigned long getID() const override { return _id; }
	/** set molecule ID */
	void setid(unsigned long id) override { _id = id; }
	/** set the molecule's component */
	void setComponent(Component *component) override { _component = component; this->updateMassInertia();}
	/** return pointer to component to which the molecule belongs */
	Component* component() const override { return _component; }
	/** get component lookUpID */
	unsigned getComponentLookUpID() const override { return _component->getLookUpId();}
	/** get position coordinate */
	double r(unsigned short d) const override { return _r[d]; }
	/** set position coordinate */
	void setr(unsigned short d, double r) override { _r[d] = r; }
	/** get velocity coordinate */
	double v(unsigned short d) const override { return _v[d]; }
	/** set velocity */
	void setv(unsigned short d, double v) override { _v[d] = v; }
	/** get molecule's mass */
	double mass() const override { return _m; }

	void setF(unsigned short d, double F) override { _F[d] = F; }
	/** get coordinate of current force onto molecule */
	double F(unsigned short d) const override {return _F[d]; }

	/** get molecule's orientation */
	const Quaternion& q() const override { return _q; }


	/** set molecule's orientation */
	void setq(Quaternion q) override{ _q = q; }

	/** get coordinate of the rotatational speed */
	double D(unsigned short d) const override { return _L[d]; }

	/** get coordinate of the current angular momentum  onto molecule */
	double M(unsigned short d) const override { return _M[d]; }

	/** get the virial **/
	double Vi(unsigned short d) const override { return _Vi[d];}
	double ViSph(unsigned short d) const override { return _ViSph[d];}
	double ViN() const override { return this->_ViSph[0];}
	double ViT() const override { return this->_ViSph[1];}


	/** get the constant correction of potential energy */
	double UpotConstCorr() const { return _upotConstCorr; }
	/** get the constant correction of one virial element */
	double ViConstCorr() const { return _ViConstCorr; }
	/** get the constant correction of NT virial element */
	double ViNConstCorr() const { return _VirNConstCorr; }
	double ViTConstCorr() const { return _VirTConstCorr; }



	void setD(unsigned short d, double D) override { this->_L[d] = D; }

	inline void move(int d, double dr) override { _r[d] += dr; }

	/** get the moment of inertia of a particle */
	double getI(unsigned short d) const override { return _I[d]; }
	/** update mass and moment of inertia by component definition */
	void updateMassInertia() override {
		if (_component != nullptr) {
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
	double v2() const override {return _v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2]; }

	/** calculate and return the square angular momentum */
	double L2() const override {return _L[0]*_L[0]+_L[1]*_L[1]+_L[2]*_L[2]; }

	/** calculate and return the square force */
	double F2() const override {return _F[0]*_F[0]+_F[1]*_F[1]+_F[2]*_F[2]; }

	/** calculate and return the square torque */
	double M2() const override {return _M[0]*_M[0]+_M[1]*_M[1]+_M[2]*_M[2]; }

	/** return the translational energy of the molecule */
	double U_trans() const override { return 0.5 * _m * v2(); }
	double U_trans_2() const override { return _m * v2(); }
	/** return the rotational energy of the molecule */
	double U_rot() override ;
	double U_rot_2() override ;
	/** return total kinetic energy of the molecule */
	double U_kin() override { return U_trans() + U_rot(); }

	void setupSoACache(CellDataSoABase * s, unsigned iLJ, unsigned iC, unsigned iD, unsigned iQ) override;

	void setSoA(CellDataSoABase * s) override;
	void setStartIndexSoA_LJ(unsigned i) override {_soa_index_lj = i;}
	void setStartIndexSoA_C(unsigned i) override {_soa_index_c = i;}
	void setStartIndexSoA_D(unsigned i) override {_soa_index_d = i;}
	void setStartIndexSoA_Q(unsigned i) override {_soa_index_q = i;}

	/* TODO: Maybe we should better do this using the component directly?
	 * In the GNU STL vector.size() causes two memory accesses and one subtraction!
	 */
	/** get number of sites */
	unsigned int numSites() const override { return _component->numSites(); }
	unsigned int numLJcenters() const override { return _component->numLJcenters(); }
	unsigned int numCharges() const override { return _component->numCharges(); }
	unsigned int numDipoles() const override { return _component->numDipoles(); }
	unsigned int numQuadrupoles() const override { return _component->numQuadrupoles(); }

	std::array<double, 3> site_d(unsigned int i) const override {
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
		} else { mardyn_assert(i < n4);
			return quadrupole_d(i - n3);
		}
	}
	std::array<double, 3> ljcenter_d(unsigned int i) const override {
		return computeLJcenter_d(i);
	}
	std::array<double, 3> charge_d(unsigned int i) const override {
		return computeCharge_d(i);
	}
	std::array<double, 3> dipole_d(unsigned int i) const override {
		return computeDipole_d(i);
	}
	std::array<double, 3> quadrupole_d(unsigned int i) const override {
		return computeQuadrupole_d(i);
	}

	std::array<double, 3> site_d_abs(unsigned int i) const override {
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
		} else { mardyn_assert(i < n4);
			return quadrupole_d_abs(i - n3);
		}
	}
	std::array<double, 3> ljcenter_d_abs(unsigned int i) const override;
	std::array<double, 3> charge_d_abs(unsigned int i) const override;
	std::array<double, 3> dipole_d_abs(unsigned int i) const override;
	std::array<double, 3> quadrupole_d_abs(unsigned int i) const override;

	std::array<double, 3> dipole_e(unsigned int i) const override;
	std::array<double, 3> quadrupole_e(unsigned int i) const override;

	std::array<double, 3> site_F(unsigned int i) const override {
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
		} else { mardyn_assert(i < n4);
			return quadrupole_F(i - n3);
		}
	}
	std::array<double, 3> ljcenter_F(unsigned int i) const override;
	std::array<double, 3> charge_F(unsigned int i) const override;
	std::array<double, 3> dipole_F(unsigned int i) const override;
	std::array<double, 3> quadrupole_F(unsigned int i) const override;

	void normalizeQuaternion() override {
		_q.normalize();
	}
	std::array<double, 3> computeLJcenter_d(unsigned int i) const override {
		mardyn_assert(_q.isNormalized());
		return _q.rotate(_component->ljcenter(i).r());
	}
	std::array<double, 3> computeCharge_d(unsigned int i) const override {
		mardyn_assert(_q.isNormalized());
		return _q.rotate(_component->charge(i).r());
	}
	std::array<double, 3> computeDipole_d(unsigned int i) const override {
		mardyn_assert(_q.isNormalized());
		return _q.rotate(_component->dipole(i).r());
	}
	std::array<double, 3> computeQuadrupole_d(unsigned int i) const override {
		mardyn_assert(_q.isNormalized());
		return _q.rotate(_component->quadrupole(i).r());
	}
	std::array<double, 3> computeDipole_e(unsigned int i) const override {
		mardyn_assert(_q.isNormalized());
		return _q.rotate(_component->dipole(i).e());
	}
	std::array<double, 3> computeQuadrupole_e(unsigned int i) const override {
		mardyn_assert(_q.isNormalized());
		return _q.rotate(_component->quadrupole(i).e());
	}


	/**
	 * get the total object memory size, together with all its members
	 * \Note You can retrieve the size of the molecule class itself simply
	 *       with the sizeof()-operator.
	 */
	unsigned long totalMemsize() const override;


	/* TODO: Is this realy necessary? Better use a function like dist2(m1.r(), m2.r()). */
	/** Calculate the difference vector and return the square (euclidean) distance.
	 *
	 *  \param molecule2 molecule to which the distance shall be calculated
	 */
	double dist2(const FullMolecule& molecule2, double dr[3]) const {
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
	void setF(double F[3]) override { for(int d = 0; d < 3; d++ ) { _F[d] = F[d]; } }

	/** set momentum acting on molecule
	 * @param M force vector (x,y,z)
	 */
	void setM(double M[3]) override { for(int d = 0; d < 3; d++ ) { _M[d] = M[d]; } }
	void setVi(double Vi[3]) override 
	{ 
		for(int d = 0; d < 3; d++) 
		{ 
			_Vi[d] = Vi[d]; 
		} 
	}
	void setViSph(double ViSph[3]) override { for(int d = 0; d < 3; d++) { _ViSph[d] = ViSph[d]; } }
	void setViN(double ViN) override { _ViSph[0] = ViN; } 
	void setViT(double ViT) override { _ViSph[1] = ViT; }

	void setUConstCorr(const double a) override { _upotConstCorr = a; }
	void setViConstCorr(const double a) override { _ViConstCorr = a/3; } // Correction term assigned to the 3 diagonal elements

	void setViNConstCorr(const double vir) { _VirNConstCorr = vir; }
	void setViTConstCorr(const double vir) { _VirTConstCorr = vir; }

	void Uadd(const double upot) override { _upot += upot; }
	void setU(const double upot) override { _upot = upot; }

	/** get the constant correction of potential energy */
	void UpotConstCorradd(double a)  { _upotConstCorr += a; }
	/** get the constant correction of one virial element */
	void ViConstCorradd(double a)  { _ViConstCorr += a; }
	/** get the constant correction of NT virial element */
	void ViNConstCorradd(double a)  { _VirNConstCorr += a; }
	void ViTConstCorradd(double a)  { _VirTConstCorr += a; }

	void Fadd(const double a[]) override { for(unsigned short d=0;d<3;++d) _F[d]+=a[d]; }
	void Madd(const double a[]) override { for(unsigned short d=0;d<3;++d) _M[d]+=a[d]; }
	void Viadd(const double a[]) override { for(unsigned short d=0;d<3;++d) _Vi[d]+=a[d]; }
	void ViSphadd(const double a[]) override { 
		for(unsigned short d=0;d<3;++d){	
		_ViSph[d]+=a[d]; }
		} 
	void ViNadd(const double a) override { 
		_ViSph[0]+=a; }
	void ViTadd(const double a) override { 
		_ViSph[1]+=a; }
	void vadd(const double ax, const double ay, const double az) override {
		_v[0] += ax; _v[1] += ay; _v[2] += az;
	}
	void vsub(const double ax, const double ay, const double az) override {
		_v[0] -= ax; _v[1] -= ay; _v[2] -= az;
	}

	void Fljcenteradd(unsigned int i, double a[]) override;
	void Fljcentersub(unsigned int i, double a[]) override;
	void Fchargeadd(unsigned int i, double a[]) override;
	void Fchargesub(unsigned int i, double a[]) override;
	void Fdipoleadd(unsigned int i, double a[]) override;
	void Fdipolesub(unsigned int i, double a[]) override;
	void Fquadrupoleadd(unsigned int i, double a[]) override;
	void Fquadrupolesub(unsigned int i, double a[]) override;

	/** First step of the leap frog integrator */
	void upd_preF(double dt) override;
	/** second step of the leap frog integrator */
	void upd_postF(double dt_halve, double& summv2, double& sumIw2) override;

	/** @brief Calculate twice the translational and rotational kinetic energies
	 * @param[out] summv2   twice the translational kinetic energy \f$ m v^2 \f$
	 * @param[out] sumIw2   twice the rotational kinetic energy \f$ I \omega^2 \f$
	 */
	void calculate_mv2_Iw2(double& summv2, double& sumIw2) override;
	void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) override;

	/**
	 * @return format of the function write(...)
	 */
	static std::string getWriteFormat();

	/** write information to stream */
	void write(std::ostream& ostrm) const override;

	/** write binary information to stream */
	void writeBinary(std::ostream& ostrm) const override;

	/** clear forces and moments */
	void clearFM() override;
	/** clear forces and moments */
	void clearVirial() override;
	/** calculate forces and moments for already given site forces */
	void calcFM() override;

	/** perform data consistency check for the molecule (only debug mode) */
	void check(unsigned long id) override;

	/** In almost all cases, molecule's caches are stored in SoAs.
	 * In some rare instances (e.g. ParticleContainer::getEnergy())
	 * a molecule should rather better exist alone and not be part of a particleCell.
	 * This function allocates a new SoA.
	 * Remember to release it when no longer necessary! */
	void buildOwnSoA() override;

	/** See above comment.*/
	void releaseOwnSoA() override;

protected:
	/** calculate forces and moments for already given site forces, for this precise site */
	void calcFM_site(const std::array<double, 3>& d, const std::array<double, 3>& F);

    Component *_component;  /**< IDentification number of its component type */
	double _r[3];  /**< position coordinates */
	double _F[3];  /**< forces */
	double _v[3];  /**< velocity */
	Quaternion _q; /**< angular orientation */
	double _M[3];  /**< torsional moment */
	double _L[3];  /**< angular momentum */
	double _Vi[3]; /** Virial tensor **/
	double _ViSph[3]; /** Spherical Virial (only has N and T component, _ViSph[2] is unused) **/
    unsigned long _id;  /**< IDentification number of that molecule */

	double _upot; /**< potential energy */
	double _upotConstCorr; /** Correction of potential energy, not changing during simulation (homogeneous system) **/
	double _ViConstCorr; /** Const. LRC of one virial element, not changing during simulation (homogeneous system) **/
	double _VirNConstCorr; /** Const. LRC of normal virial element **/
	double _VirTConstCorr; /** Const. LRC of normal virial element **/




	double _m; /**< total mass */
	double _I[3]{0.,0.,0.},_invI[3]{0.,0.,0.};  // moment of inertia for principal axes and it's inverse

	/* absolute positions are stored in the soa. Work-arounds to get the relative ones*/
	CellDataSoA * _soa;
	unsigned _soa_index_lj;
	unsigned _soa_index_c;
	unsigned _soa_index_d;
	unsigned _soa_index_q;

};

std::ostream& operator<<( std::ostream& os, const FullMolecule& m );

#endif /* FULLMOLECULE_H_ */
