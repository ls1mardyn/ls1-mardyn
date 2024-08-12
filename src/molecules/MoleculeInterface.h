/*
 * MoleculeInterface.h
 *
 *  Created on: 21 Jan 2017
 *      Author: tchipevn
 */

#ifndef SRC_MOLECULES_MOLECULEINTERFACE_H_
#define SRC_MOLECULES_MOLECULEINTERFACE_H_

#include <array>

#include "molecules/Component.h"
#include "molecules/Quaternion.h"
#include "particleContainer/adapter/CellDataSoABase.h"
#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"

class MoleculeInterface {
public:
	static std::array<vcp_real_calc, 3> convert_double_to_vcp_real_calc(const std::array<double,3>& v) {
		std::array<vcp_real_calc, 3> ret;
		for (int d = 0; d < 3; ++d) {
			ret[d] = static_cast<vcp_real_calc>(v[d]);
		}
		return ret;
	}

	virtual ~MoleculeInterface() {}
	virtual unsigned long getID() const = 0;
	virtual void setid(unsigned long id) = 0;
	virtual void setComponent(Component *component) = 0;
	virtual void setr(unsigned short d, double r) = 0;
	virtual void setv(unsigned short d, double v) = 0;
	virtual void setF(unsigned short d, double F) = 0;
	unsigned int componentid() const {
		return component()->ID();
	}
	virtual Component* component() const = 0;
	virtual unsigned getComponentLookUpID() const {
		return component()->getLookUpId();
	}
	virtual double r(unsigned short d) const = 0;
	std::array<double, 3> r_arr() const {
		std::array<double, 3> ret;
		for(int d=0; d < 3; ++d) {
			ret[d] = r(d);
		}
		return ret;
	}
	virtual double v(unsigned short d) const = 0;
	std::array<double, 3> v_arr() const {
		std::array<double, 3> ret = {v(0), v(1), v(2)};
		return ret;
	}
	virtual double mass() const {
		return component()->m();
	}
	virtual double F(unsigned short d) const = 0;

	virtual std::array<double, 3> F_arr() {
		std::array<double, 3> f_ret{F(0), F(1), F(2)};
		return f_ret;
	}

	virtual std::array<double, 3> M_arr() {
		std::array<double, 3> m_ret{M(0), M(1), M(2)};
		return m_ret;
	}
	virtual std::array<double, 3> Vi_arr() {
		std::array<double, 3> vi_ret{Vi(0), Vi(1), Vi(2)};
		return vi_ret;
	}
	virtual std::array<double, 3> ViSph_arr() {
		std::array<double, 3> viSph_ret{ViSph(0), ViSph(1), ViSph(2)};
		return viSph_ret;
	}

	virtual const Quaternion& q() const = 0;


	virtual void setq(Quaternion q)= 0;

	virtual double D(unsigned short d) const = 0;
	std::array<double, 3> D_arr() const {
		std::array<double, 3> ret{D(0), D(1), D(2)};
		return ret;
	}
	virtual double M(unsigned short d) const = 0;
	virtual double Vi(unsigned short d) const = 0;
	virtual double ViSph(unsigned short d) const = 0;
	virtual double ViN() const = 0;
	virtual double ViT() const = 0;


	// /** get the constant correction of potential energy */
	// virtual double UpotConstCorr() const = 0;
	// /** get the constant correction of one virial element */
	// virtual double ViConstCorr() const = 0;
	// /** get the constant correction of NT virial element */
	// virtual double ViNConstCorr() const = 0;
	// virtual double ViTConstCorr() const = 0;

	virtual void setD(unsigned short d, double D) = 0;

	inline virtual void move(int d, double dr) = 0;

	//by Stefan Becker
	virtual double getI(unsigned short d) const = 0;


	virtual double v2() const {return v(0)*v(0)+v(1)*v(1)+v(2)*v(2); }
	virtual double F2() const {return F(0)*F(0)+F(1)*F(1)+F(2)*F(2); }
	virtual double L2() const {return D(0)*D(0)+D(1)*D(1)+D(2)*D(2); }
	virtual double M2() const {return M(0)*M(0)+M(1)*M(1)+M(2)*M(2); }

	virtual double U_trans() const { return 0.5 * mass() * v2(); }
	virtual double U_trans_2() const { return mass() * v2(); }
	virtual double U_rot() = 0;
	virtual double U_rot_2() = 0;
	virtual double U_kin() { return U_trans() + U_rot(); }

	virtual void updateMassInertia() = 0;

	virtual void setupSoACache(CellDataSoABase * const s, unsigned iLJ, unsigned iC, unsigned iD, unsigned iQ) = 0;

	virtual void setSoA(CellDataSoABase * const s) = 0;
	virtual void setStartIndexSoA_LJ(unsigned i) = 0;
	virtual void setStartIndexSoA_C(unsigned i) = 0;
	virtual void setStartIndexSoA_D(unsigned i) = 0;
	virtual void setStartIndexSoA_Q(unsigned i) = 0;

	virtual unsigned int numSites() const = 0;
	virtual unsigned int numLJcenters() const = 0;
	virtual unsigned int numCharges() const = 0;
	virtual unsigned int numDipoles() const = 0;
	virtual unsigned int numQuadrupoles() const = 0;

	// relative positions to CM
	virtual std::array<double, 3> site_d(unsigned int i) const = 0;
	virtual std::array<double, 3> ljcenter_d(unsigned int i) const = 0;
	virtual std::array<double, 3> charge_d(unsigned int i) const = 0;
	virtual std::array<double, 3> dipole_d(unsigned int i) const = 0;
	virtual std::array<double, 3> quadrupole_d(unsigned int i) const = 0;

	// absolute positions to global zero
	virtual std::array<double, 3> site_d_abs(unsigned int i) const = 0;
	virtual std::array<double, 3> ljcenter_d_abs(unsigned int i) const = 0;
	virtual std::array<double, 3> charge_d_abs(unsigned int i) const = 0;
	virtual std::array<double, 3> dipole_d_abs(unsigned int i) const = 0;
	virtual std::array<double, 3> quadrupole_d_abs(unsigned int i) const = 0;

	virtual std::array<double, 3> dipole_e(unsigned int i) const = 0;
	virtual std::array<double, 3> quadrupole_e(unsigned int i) const = 0;

	virtual std::array<double, 3> site_F(unsigned int i) const = 0;
	virtual std::array<double, 3> ljcenter_F(unsigned int i) const = 0;
	virtual std::array<double, 3> charge_F(unsigned int i) const = 0;
	virtual std::array<double, 3> dipole_F(unsigned int i) const = 0;
	virtual std::array<double, 3> quadrupole_F(unsigned int i) const = 0;

	virtual void normalizeQuaternion() = 0;
	virtual std::array<double, 3> computeLJcenter_d(unsigned int i) const = 0;
	virtual std::array<double, 3> computeCharge_d(unsigned int i) const = 0;
	virtual std::array<double, 3> computeDipole_d(unsigned int i) const = 0;
	virtual std::array<double, 3> computeQuadrupole_d(unsigned int i) const = 0;
	virtual std::array<double, 3> computeDipole_e(unsigned int i) const = 0;
	virtual std::array<double, 3> computeQuadrupole_e(unsigned int i) const = 0;


	virtual unsigned long totalMemsize() const = 0;

	double dist2(const MoleculeInterface& molecule2, double dr[3]) const {
		double d2 = 0.;
		for (unsigned short d = 0; d < 3; d++) {
			dr[d] = molecule2.r(d) - r(d);
			d2 += dr[d] * dr[d];
		}
		return d2;
	}

	//calculates orientation angle for ARDF
	double orientationAngle(const MoleculeInterface& molecule2, double dr[3], double d2) const {

		double cosPhi = 0.;
		double orientationVector[3];
		double orientationVectorSquared = 0.;
		double roundingThreshold = 0.0001;

		orientationVector[0] = 2. * (q().qx() * q().qz() + q().qw() * q().qy());
		orientationVector[1] = 2. * (q().qy() * q().qz() - q().qw() * q().qx());
		orientationVector[2] = 1. - 2. * (q().qx() * q().qx() + q().qy() * q().qy());



		for (unsigned short d = 0; d < 3; d++) {
			dr[d] = molecule2.r(d) - r(d);
			orientationVectorSquared += orientationVector[d] * orientationVector[d];
		}

		for (unsigned short d = 0; d < 3; d++) {
			cosPhi += orientationVector[d] * dr[d] / sqrt(orientationVectorSquared) / sqrt(d2);
		}
		return cosPhi;

	}

	virtual void setF(double F[3]) = 0;
	virtual void setM(double M[3]) = 0;
	virtual void setVi(double Vi[3]) = 0;
	virtual void setViSph(double Vi[3]) = 0;
	virtual void setViN(double ViN) = 0;
	virtual void setViT(double ViT) = 0;
	void scale_v(double s) {
		for(int d = 0; d < 3; ++d) {
			setv(d, v(d) * s);
		}
	}
	void scale_v(double s, double offx, double offy, double offz) {
		vsub(offx, offy, offz);
		scale_v(s);
		vadd(offx, offy, offz);
	}
	void scale_F(double s) {
		double Fscaled[3] = {F(0) * s, F(1) * s, F(2) * s};
		setF(Fscaled);
	}
	void scale_D(double s) {
		for (int d = 0; d < 3; ++d) {
			setD(d, D(d) * s);
		}
	}
	void scale_M(double s) {
		double Mscaled[3] = {M(0) * s, M(1) * s, M(2) * s};
		setM(Mscaled);
	}
	virtual void Fadd(const double a[]) = 0;
	virtual void Madd(const double a[]) = 0;
	virtual void Viadd(const double a[]) = 0;
	virtual void ViSphadd(const double a[]) = 0;
	virtual void ViNadd(const double a) = 0;
	virtual void ViTadd(const double a) = 0;
	virtual void vadd(const double ax, const double ay, const double az) = 0;
	virtual void vsub(const double ax, const double ay, const double az) = 0;

	virtual void setUConstCorr(const double a) = 0;
	virtual void setViConstCorr(const double a) = 0;

	virtual void Uadd(const double upot) = 0;
	virtual void setU(const double upot) = 0;
	
	virtual void Fljcenteradd(unsigned int i, double a[]) = 0;
	virtual void Fljcentersub(unsigned int i, double a[]) = 0;
	virtual void Fchargeadd(unsigned int i, double a[]) = 0;
	virtual void Fchargesub(unsigned int i, double a[]) = 0;
	virtual void Fdipoleadd(unsigned int i, double a[]) = 0;
	virtual void Fdipolesub(unsigned int i, double a[]) = 0;
	virtual void Fquadrupoleadd(unsigned int i, double a[]) = 0;
	virtual void Fquadrupolesub(unsigned int i, double a[]) = 0;

	// Leapfrog integration:
	virtual void upd_preF(double dt) = 0;
	virtual void upd_postF(double dt_halve, double& summv2, double& sumIw2) = 0;

	// Integration for single-centered molecules in RMM mode.
	// Explicit Euler and LeapFrog without time-splitting for single-centered molecules are identical
	// up to interpretation of whether the velocities are stored at (t) or at (t+dt/2)  (right?)
	// Placing it here, so that we can
	// verify the RMM code against the non-RMM code.
	void ee_upd_preF(double dt) {
		for (unsigned short d = 0; d < 3; ++d) {
			setr(d,r(d) + dt * v(d));
		}
	}
	void ee_upd_postF(double
#ifndef ENABLE_REDUCED_MEMORY_MODE
		dt
#endif
		, double& summv2) {

		//calcFM(); //NOTE: This was moved to simulation.cpp calculateForces() and is called in simulate()

#ifndef ENABLE_REDUCED_MEMORY_MODE
		double dtInvM = dt / component()->m();
		for (unsigned short d = 0; d < 3; ++d) {
			setv(d, v(d) + dtInvM * F(d));
		}
#endif /*ENABLE_REDUCED_MEMORY_MODE */

		summv2 += component()->m() * v2();
	}

	virtual void calculate_mv2_Iw2(double& summv2, double& sumIw2) = 0;
	virtual void calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) = 0;
	static std::string getWriteFormat(); // TODO
	virtual void write(std::ostream& ostrm) const = 0;
	virtual void writeBinary(std::ostream& ostrm) const = 0;
	virtual void clearFM() = 0;
	virtual void clearVirial() = 0;
	virtual void calcFM() = 0;
	virtual void check(unsigned long id) = 0;

	//! NOTE: Vectorized force calculation doesn't use this anymore
	//!
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
	bool isLessThan(const MoleculeInterface& m2) const;

	/**
	 * \brief test whether molecule is inside a cuboid region
	 * @param l lower left front corner of cube (equality allowed)
	 * @param u upper right back corner of cube (equality not allowed)
	 * @return true if molecule is contained in the box, false otherwise
	 */
	virtual bool inBox(const double l[3], const double u[3]) const {
		bool in = true;
		for (int d = 0; d < 3; ++d) {
#ifdef __INTEL_COMPILER
			#pragma float_control(precise, on)
			#pragma fenv_access(on)
#endif
			in &= (r(d) >= l[d] and r(d) < u[d]);
		}
		return in;
	}

	/** In almost all cases, molecule's caches are stored in SoAs.
	 * In some rare instances (e.g. ParticleContainer::getEnergy())
	 * a molecule should rather better exist alone and not be part of a particleCell.
	 * This function allocates a new SoA.
	 * Remember to release it when no longer necessary! */
	virtual void buildOwnSoA() = 0;

	/** See above comment.*/
	virtual void releaseOwnSoA() = 0;
};

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

#endif /* SRC_MOLECULES_MOLECULEINTERFACE_H_ */
