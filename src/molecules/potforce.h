#ifndef POTFORCE_H_
#define POTFORCE_H_

/**  @file  */

#include "molecules/Comp2Param.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "Domain.h"

/**
 * The formulas for dipole-dipole, dipole-quadrupole and quadrupole-quadrupole are from
 * Gray and Gubbins, Theory of molecular Fluids, Oxford, 1984.
 * The individual formulas are specified below.
 * They are from Section "Explicit angle dependence" in Chapter 2.4 Multipole Interactions.
 * They are given in electrostatic units, meaning the 1 / 4 pi eps_0 Coulomb constant is left out.
 */

/** @brief Calculates potential and force between 2 Lennard-Jones 12-6 centers. */
inline void PotForceLJ(const double dr[3], const double& dr2,
                       const double& eps24, const double& sig2,
                       double f[3], double& u6)
{
	double invdr2 = 1. / dr2;
	double lj6 = sig2 * invdr2; lj6 = lj6 * lj6 * lj6;
	double lj12 = lj6 * lj6;
	double lj12m6 = lj12 - lj6;
	u6 = eps24 * lj12m6;
	double fac = eps24 * (lj12 + lj12m6) * invdr2;
	for (unsigned short d = 0; d < 3; ++d)
		f[d] = fac * dr[d];
}


/** @brief Calculate potential and force between 2 Dipoles.
 * Formula (2.180) in Gray and Gubbins, see comment at start of this file.
 */
inline void PotForce2Dipole(const double dr[3], const double& dr2, const double* eii, const double* ejj,
                            const double& my2, const double& rffac,
                            double f[3], double m1[3], double m2[3], double& u, double& MyRF)
{
	double invdr2 = 1. / dr2;
	double invdr1 = sqrt(invdr2);
	double myfac = my2 * invdr2 * invdr1;
	double costi = 0.;
	double costj = 0.;
	double cosgij = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		const double& drd = dr[d];
		const double& eiid = eii[d];
		const double& ejjd = ejj[d];
		costi += eiid * drd;
		costj += ejjd * drd;
		cosgij += eiid * ejjd;
	}
	costi *= invdr1;
	costj *= invdr1;
	u = myfac * (cosgij - 3. * costi * costj);
	MyRF -= rffac * cosgij;
	const double partialRijInvdr1 = -3. * u * invdr2;
	const double partialTiInvdr1 = -myfac * 3. * costj * invdr1;
	const double partialTjInvdr1 = -myfac * 3. * costi * invdr1;
	const double& partialGij = myfac;
	const double fac = -partialRijInvdr1 + (costi*partialTiInvdr1 + costj*partialTjInvdr1) * invdr1;
	for (unsigned short d = 0; d < 3; ++d)
		f[d] = fac * dr[d] - partialTiInvdr1*eii[d] - partialTjInvdr1*ejj[d];

	double eiXej[3], eXrij[3];
	eiXej[0] = eii[1] * ejj[2] - eii[2] * ejj[1];
	eiXej[1] = eii[2] * ejj[0] - eii[0] * ejj[2];
	eiXej[2] = eii[0] * ejj[1] - eii[1] * ejj[0];
	eXrij[0] = eii[1] * dr[2] - eii[2] * dr[1];
	eXrij[1] = eii[2] * dr[0] - eii[0] * dr[2];
	eXrij[2] = eii[0] * dr[1] - eii[1] * dr[0];
	for (unsigned short d = 0; d < 3; ++d)
		m1[d] = -partialTiInvdr1 * eXrij[d] + (-partialGij + rffac) * eiXej[d];
	eXrij[0] = ejj[1] * dr[2] - ejj[2] * dr[1];
	eXrij[1] = ejj[2] * dr[0] - ejj[0] * dr[2];
	eXrij[2] = ejj[0] * dr[1] - ejj[1] * dr[0];
	for (unsigned short d = 0; d < 3; ++d)
		m2[d] = -partialTjInvdr1 * eXrij[d] + (partialGij - rffac) * eiXej[d];
}


/** @brief Calculate potential and force between 2 Quadrupoles.
 * Formula (2.184) in Gray and Gubbins, see comment at start of this file.
 */
inline void PotForce2Quadrupole(const double dr[3], const double& dr2, const double* eii, const double* ejj,
                                const double& q2075,
                                double f[3], double m1[3], double m2[3], double& u)
{
	double invdr2 = 1. / dr2;
	double invdr1 = sqrt(invdr2);
	double qfac = q2075 * invdr2 * invdr2 * invdr1;
	double costi = 0.;
	double costj = 0.;
	double cosgij = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		const double& drd = dr[d];
		const double& eiid = eii[d];
		const double& ejjd = ejj[d];
		costi += eiid * drd;
		costj += ejjd * drd;
		cosgij += eiid * ejjd;
	}
	costi *= invdr1;
	costj *= invdr1;
	double cos2ti = costi * costi;
	double cos2tj = costj * costj;
	double term = (cosgij - 5. * costi * costj);
	u = qfac * (1. - 5. * (cos2ti + cos2tj) - 15. * cos2ti * cos2tj + 2. * term * term);
	const double partialRijInvdr1 = -5. * u * invdr2;
	const double partialTiInvdr1 = -qfac * 10. * (costi + 3. * costi * cos2tj + 2. * costj * term) * invdr1;
	const double partialTjInvdr1 = -qfac * 10. * (costj + 3. * cos2ti * costj + 2. * costi * term) * invdr1;
	const double partialGij = qfac * 4. * term;
	const double fac = -partialRijInvdr1 + (costi*partialTiInvdr1 + costj*partialTjInvdr1) * invdr1;
	for (unsigned short d = 0; d < 3; ++d)
		f[d] = fac * dr[d] - partialTiInvdr1*eii[d] - partialTjInvdr1*ejj[d];

	double eiXej[3], eXrij[3];
	eiXej[0] = eii[1] * ejj[2] - eii[2] * ejj[1];
	eiXej[1] = eii[2] * ejj[0] - eii[0] * ejj[2];
	eiXej[2] = eii[0] * ejj[1] - eii[1] * ejj[0];
	eXrij[0] = eii[1] * dr[2] - eii[2] * dr[1];
	eXrij[1] = eii[2] * dr[0] - eii[0] * dr[2];
	eXrij[2] = eii[0] * dr[1] - eii[1] * dr[0];
	for (unsigned short d = 0; d < 3; ++d)
		m1[d] = -partialTiInvdr1 * eXrij[d] - partialGij * eiXej[d];

	eXrij[0] = ejj[1] * dr[2] - ejj[2] * dr[1];
	eXrij[1] = ejj[2] * dr[0] - ejj[0] * dr[2];
	eXrij[2] = ejj[0] * dr[1] - ejj[1] * dr[0];
	for (unsigned short d = 0; d < 3; ++d)
		m2[d] = -partialTjInvdr1 * eXrij[d] + partialGij * eiXej[d];
}


/** @brief Calculate potential and force between a Dipole and Quadrupole.
 *  Formula (2.182) in Gray and Gubbins, see comment at start of this file.
 */
inline void PotForceDiQuadrupole(const double dr[3], const double& dr2, const double* eii, const double* ejj,
                                 const double& myq15,
                                 double f[3], double m1[3], double m2[3], double& u)
{
	double invdr2 = 1. / dr2;
	double invdr1 = sqrt(invdr2);
	double myqfac = myq15 * invdr2 * invdr2;
	double costi = 0.;
	double costj = 0.;
	double cosgij = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		const double& drd = dr[d];
		const double& eiid = eii[d];
		const double& ejjd = ejj[d];
		costi += eiid * drd;
		costj += ejjd * drd;
		cosgij += eiid * ejjd;
	}
	costi *= invdr1;
	costj *= invdr1;
	double cos2tj = costj * costj;
	u = myqfac * (-costi * (5. * cos2tj - 1.) + 2. * cosgij * costj);
	const double partialRijInvdr1 = -4. * u * invdr2;
	const double partialTiInvdr1 = myqfac * (-5. * cos2tj + 1.) * invdr1;
	const double partialTjInvdr1 = myqfac * 2. * (-5. * costi * costj + cosgij) * invdr1;
	const double partialGij = myqfac * 2. * costj;
	const double fac = -partialRijInvdr1 + (costi*partialTiInvdr1 + costj*partialTjInvdr1) * invdr1;
	for (unsigned short d = 0; d < 3; ++d)
		f[d] = fac * dr[d] - partialTiInvdr1*eii[d] - partialTjInvdr1*ejj[d];

	double eiXej[3], eXrij[3];
	eiXej[0] = eii[1] * ejj[2] - eii[2] * ejj[1];
	eiXej[1] = eii[2] * ejj[0] - eii[0] * ejj[2];
	eiXej[2] = eii[0] * ejj[1] - eii[1] * ejj[0];
	eXrij[0] = eii[1] * dr[2] - eii[2] * dr[1];
	eXrij[1] = eii[2] * dr[0] - eii[0] * dr[2];
	eXrij[2] = eii[0] * dr[1] - eii[1] * dr[0];
	for (unsigned short d = 0; d < 3; ++d)
		m1[d] = -partialTiInvdr1 * eXrij[d] - partialGij * eiXej[d];

	eXrij[0] = ejj[1] * dr[2] - ejj[2] * dr[1];
	eXrij[1] = ejj[2] * dr[0] - ejj[0] * dr[2];
	eXrij[2] = ejj[0] * dr[1] - ejj[1] * dr[0];
	for (unsigned short d = 0; d < 3; ++d)
		m2[d] = -partialTjInvdr1 * eXrij[d] + partialGij * eiXej[d];
}


/** @brief Calculate potential and force between two point charges. 
 * Coulomb's law.
 */
inline void PotForce2Charge(const double dr[3], const double& dr2,
                            const double& q1q2per4pie0, double f[3], double& u)
{
	double invdr2 = 1.0 / dr2;
	double invdr = sqrt(invdr2);
	u = q1q2per4pie0 * invdr;
	const double fac = u * invdr2;
	for (unsigned short d = 0; d < 3; d++)
		f[d] = fac * dr[d];
}


/** @brief Calculate potential and force between an electric charge and a quadrupole.
 * Just Quadrupole electric field times charge magnitude.
 */
inline void PotForceChargeQuadrupole(const double dr[3], const double& dr2,
                                     const double* ejj, const double& qQ05per4pie0,
                                     double f[3], double m2[3], double& u)
{
	double invdr2 = 1.0 / dr2;
	double invdr = sqrt(invdr2);
	double costj = 0;
	for (unsigned short d = 0; d < 3; d++) {
		costj += ejj[d] * dr[d];
	}
	costj *= invdr;
	double qQinv4dr3 = qQ05per4pie0 * invdr * invdr2;
	u = qQinv4dr3 * (3.0 * costj * costj - 1);

	const double partialRijInvdr1 = -3.0 * u * invdr2;
	const double partialTjInvdr1 = 6.0 * costj * qQinv4dr3 * invdr;
	const double fac = costj * partialTjInvdr1 * invdr - partialRijInvdr1;
	for (unsigned short d = 0; d < 3; ++d)
		f[d] = fac*dr[d] - partialTjInvdr1 * ejj[d];

	double minuseXrij[3];
	minuseXrij[0] = ejj[2]*dr[1] - ejj[1]*dr[2];
	minuseXrij[1] = ejj[0]*dr[2] - ejj[2]*dr[0];
	minuseXrij[2] = ejj[1]*dr[0] - ejj[0]*dr[1];
	for (unsigned short d = 0; d < 3; ++d)
		m2[d] = partialTjInvdr1 * minuseXrij[d];
}


/** @brief Calculate potential and force between an electric charge and a dipole.
 * Just Dipole electric field times charge magnitude.
 */
inline void PotForceChargeDipole(const double dr[3], const double& dr2,
                                 const double* ejj, const double& minusqmyper4pie0,
                                 double f[3], double m2[3], double& u)
{
	double invdr2 = 1.0 / dr2;
	double invdr = sqrt(invdr2);
	double costj = 0;
	for (unsigned short d = 0; d < 3; d++) {
		costj += ejj[d] * dr[d];
	}
	costj *= invdr;
	double uInvcostj = minusqmyper4pie0 * invdr2;
	u = uInvcostj * costj;

	// const double partialRijInvdr1 = -2.0 * u * invdr2;
	const double partialTjInvdr1 = uInvcostj*invdr;
	const double fac = 3.0 * u * invdr2;
	for (unsigned short d = 0; d < 3; ++d)
		f[d] = fac*dr[d] - partialTjInvdr1 * ejj[d];

	double minuseXrij[3];
	minuseXrij[0] = ejj[2]*dr[1] - ejj[1]*dr[2];
	minuseXrij[1] = ejj[0]*dr[2] - ejj[2]*dr[0];
	minuseXrij[2] = ejj[1]*dr[0] - ejj[0]*dr[1];
	for (unsigned short d = 0; d < 3; ++d)
		m2[d] = partialTjInvdr1 * minuseXrij[d];
}


/** @brief Calculate potential and force between two molecules including all site-site interactions.
 *
 * Calculates the potential energy and force between two molecules i and j.
 * The interaction paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
 *
 * @param[in]  mi   molecule i
 * @param[in]  mj   molecule j
 * @param[in]  params    reference to the corresponding interaction parameter stream initialized via Comp2Param::initialize
 * @param[in]  drm   distance vector from molecule j to molecule i
 * @param[out] Upot6LJ   potential energy resulting from Lennard Jones interactions
 * @param[out] UpotXpoles   potential energy resulting from Charge, Dipole and Quadrupole interactions
 * @param[out] MyRF
 * @todo Document parameter (is this reaction field?)
 * @param[out] Virial   Virial
 * @param[in]  caculateLJ    enable or disable calculation of Lennard Jones interactions
 */
inline void PotForce(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3], bool calculateLJ)
// ???better calc Virial, when molecule forces are calculated:
//    summing up molecule virials instead of site virials???
{ // Force Calculation
	double f[3];
	double u;
	double drs[3], dr2; // site distance vector & length^2
	Virial[0]=0.;
	Virial[1]=0.;
	Virial[2]=0.;
	double VNi = 0.;
	double VTi = 0.;
	double VNj = 0.;
	double VTj = 0.;

	// LJ centers
	// no LJ interaction between solid atoms of the same component

	const unsigned int nc1 = mi.numLJcenters();
	const unsigned int nc2 = mj.numLJcenters();
	for (unsigned int si = 0; si < nc1; ++si) {
		const std::array<double,3> dii = mi.ljcenter_d_abs(si);
		for (unsigned int sj = 0; sj < nc2; ++sj) {
			const std::array<double,3> djj = mj.ljcenter_d_abs(sj);
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			double eps24;
			params >> eps24;
			double sig2;
			params >> sig2;
			double shift6;
			params >> shift6; // must be 0.0 for full LJ
			if (calculateLJ) {
				PotForceLJ(drs, dr2, eps24, sig2, f, u);
				u += shift6;
				mi.Fljcenteradd(si, f);
				mj.Fljcentersub(sj, f);
				Upot6LJ += u;
				for (unsigned short d = 0; d < 3; ++d)
					Virial[d] += 0.5*drm[d] * f[d];

				// pN, pT
        double center = 0.5*global_simulation->getDomain()->getGlobalLength(0);

				double ksi_i[3];
				double ksi_j[3];
				for (unsigned short d = 0; d < 3; ++d) {
					ksi_i[d] = center - mi.r(d);
					ksi_j[d] = center - mj.r(d);
				}

				double absi2 = pow(ksi_i[0],2)+pow(ksi_i[1],2)+pow(ksi_i[2],2);
				double absj2 = pow(ksi_j[0],2)+pow(ksi_j[1],2)+pow(ksi_j[2],2);

				double drm_ni = 0.;
				double drm_nj = 0.;
				double drs_ni = 0.;
				double drs_nj = 0.;
				for (unsigned short d = 0; d < 3; ++d) {
					drm_ni += drm[d] * ksi_i[d];
					drm_nj += drm[d] * ksi_j[d];
					drs_ni += drs[d] * ksi_i[d];
					drs_nj += drs[d] * ksi_j[d];

				}
				double fac = (f[0]+f[1]+f[2])/(drs[0]+drs[1]+drs[2]);

				VNi += 0.5*fac*drm_ni*drs_ni/absi2;
				VNj += 0.5*fac*drm_nj*drs_nj/absj2;

				double drm_ti[3];
				double drs_ti[3];
				double drm_tj[3];
				double drs_tj[3];
				for (unsigned short d = 0; d < 3; ++d) {
					drm_ti[d] = drm[d] - drm_ni*ksi_i[d]/absi2;
					drs_ti[d] = drs[d] - drs_ni*ksi_i[d]/absi2;
					drm_tj[d] = drm[d] - drm_nj*ksi_j[d]/absj2;
					drs_tj[d] = drs[d] - drs_nj*ksi_j[d]/absj2;
				}
				for (unsigned short d = 0; d < 3; ++d) {
					VTi += 0.25*fac*drm_ti[d]*drs_ti[d];
					VTj += 0.25*fac*drm_tj[d]*drs_tj[d];
				}
				
			}
		}
	}

	double m1[3], m2[3]; // angular momenta

	const unsigned ne1 = mi.numCharges();
	const unsigned ne2 = mj.numCharges();
	const unsigned int nq1 = mi.numQuadrupoles();
	const unsigned int nq2 = mj.numQuadrupoles();
	const unsigned int nd1 = mi.numDipoles();
	const unsigned int nd2 = mj.numDipoles();
	for (unsigned si = 0; si < ne1; si++) {
		const std::array<double,3> dii = mi.charge_d_abs(si);
		// Charge-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = mj.charge_d_abs(sj);
			double q1q2per4pie0; // 4pie0 = 1 in reduced units
			params >> q1q2per4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);

			mi.Fchargeadd(si, f);
			mj.Fchargesub(sj, f);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial[d] += 0.5*drm[d] * f[d];
		}
		// Charge-Quadrupole
		for (unsigned sj = 0; sj < nq2; sj++) {
			const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.quadrupole_e(sj);
			PotForceChargeQuadrupole(drs, dr2, ejj.data(), qQ05per4pie0, f, m2, u);

			mi.Fchargeadd(si, f);
			mj.Fquadrupolesub(sj, f);
			mj.Madd(m2);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial[d] += 0.5*drm[d] * f[d];
		}
		// Charge-Dipole
		for (unsigned sj = 0; sj < nd2; sj++) {
			const std::array<double,3> djj = mj.dipole_d_abs(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.dipole_e(sj);
			PotForceChargeDipole(drs, dr2, ejj.data(), minusqmyper4pie0, f, m2, u);

			mi.Fchargeadd(si, f);
			mj.Fdipolesub(sj, f);
			mj.Madd(m2);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial[d] += 0.5*drm[d] * f[d];
		}
	}
	for (unsigned int si = 0; si < nq1; ++si) {
		const std::array<double,3> dii = mi.quadrupole_d_abs(si);
		const std::array<double,3> eii = mi.quadrupole_e(si);

		// Quadrupole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = mj.charge_d_abs(sj);
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForceChargeQuadrupole(drs, dr2, eii.data(), qQ05per4pie0, f, m1, u);

			mi.Fquadrupolesub(si, f);
			mj.Fchargeadd(sj, f);
			mi.Madd(m1);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial[d] -= 0.5*drm[d] * f[d];
		}
		// Quadrupole-Quadrupole -------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			//double drs[3];
			const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
			double q2075;
			params >> q2075;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.quadrupole_e(sj);
			PotForce2Quadrupole(drs, dr2, eii.data(), ejj.data(), q2075, f, m1, m2, u);

			mi.Fquadrupoleadd(si, f);
			mj.Fquadrupolesub(sj, f);
			mi.Madd(m1);
			mj.Madd(m2);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				Virial[d] += 0.5*drm[d] * f[d];
		}
		// Quadrupole-Dipole -----------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			//double drs[3];
			const std::array<double,3> djj = mj.dipole_d_abs(sj);
			double qmy15;
			params >> qmy15;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.dipole_e(sj);
			PotForceDiQuadrupole(drs, dr2, ejj.data(), eii.data(), qmy15, f, m2, m1, u);

			mi.Fquadrupolesub(si, f);
			mj.Fdipoleadd(sj, f);
			mi.Madd(m1);
			mj.Madd(m2);
			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial[d] -= 0.5*drm[d] * f[d];
		}
	}
	for (unsigned int si = 0; si < nd1; ++si) {
		const std::array<double,3> dii = mi.dipole_d_abs(si);
		const std::array<double,3> eii = mi.dipole_e(si);
		// Dipole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = mj.charge_d_abs(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForceChargeDipole(drs, dr2, eii.data(), minusqmyper4pie0, f, m1, u);

			mi.Fdipolesub(si, f);
			mj.Fchargeadd(sj, f);
			mi.Madd(m1);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial[d] -= 0.5*drm[d] * f[d];
		}
		// Dipole-Quadrupole -----------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			//double drs[3];
			const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
			double myq15;
			params >> myq15;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.quadrupole_e(sj);
			PotForceDiQuadrupole(drs, dr2, eii.data(), ejj.data(), myq15, f, m1, m2, u);

			mi.Fdipoleadd(si, f);
			mj.Fquadrupolesub(sj, f);
			mi.Madd(m1);
			mj.Madd(m2);
			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				Virial[d] += 0.5*drm[d] * f[d];
		}
		// Dipole-Dipole ---------------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			const std::array<double,3> djj = mj.dipole_d_abs(sj);
			double my2;
			params >> my2;
			double rffac;
			params >> rffac;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.dipole_e(sj);
			PotForce2Dipole(drs, dr2, eii.data(), ejj.data(), my2, rffac, f, m1, m2, u, MyRF);

			mi.Fdipoleadd(si, f);
			mj.Fdipolesub(sj, f);
			mi.Madd(m1);
			mj.Madd(m2);
			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				Virial[d] += 0.5*drm[d] * f[d];
		}
	}

	mi.Viadd(Virial);
	mj.Viadd(Virial);
	mi.VirNadd(VNi);
	mj.VirNadd(VNj);
	mi.VirTadd(VTi);
	mj.VirTadd(VTj);

	// check whether all parameters were used
	mardyn_assert(params.eos());
}

/** @brief Calculates the LJ and electrostatic potential energy of the mi-mj interaction (no multi-body potentials are considered) */
inline void FluidPot(Molecule& mi, Molecule& mj, ParaStrm& params, double /*drm*/[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, bool calculateLJ)
{
	double f[3];
	double u;
	double drs[3], dr2; // site distance vector & length^2
	// no LJ interaction between equal solid atoms

	const unsigned int nc1 = mi.numLJcenters();
	const unsigned int nc2 = mj.numLJcenters();
	for (unsigned int si = 0; si < nc1; ++si) {
		const std::array<double,3> dii = mi.ljcenter_d_abs(si);
		for (unsigned int sj = 0; sj < nc2; ++sj) {
			const std::array<double,3> djj = mj.ljcenter_d_abs(sj);
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			double eps24;
			params >> eps24;
			double sig2;
			params >> sig2;
			double shift6;
			params >> shift6; // must be 0.0 for full LJ

			if (calculateLJ) {
				PotForceLJ(drs, dr2, eps24, sig2, f, u);
				u += shift6;
				Upot6LJ += u;
			}
		}
	}


	double m1[3], m2[3]; // angular momenta

	const unsigned ne1 = mi.numCharges();
	const unsigned ne2 = mj.numCharges();
	const unsigned int nq1 = mi.numQuadrupoles();
	const unsigned int nq2 = mj.numQuadrupoles();
	const unsigned int nd1 = mi.numDipoles();
	const unsigned int nd2 = mj.numDipoles();
	for (unsigned si = 0; si < ne1; si++) {
		const std::array<double,3> dii = mi.charge_d_abs(si);
		// Charge-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = mj.charge_d_abs(sj);
			double q1q2per4pie0; // 4pie0 = 1 in reduced units
			params >> q1q2per4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);
			UpotXpoles += u;
		}
		// Charge-Quadrupole
		for (unsigned sj = 0; sj < nq2; sj++) {
			const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.quadrupole_e(sj);
			PotForceChargeQuadrupole(drs, dr2, ejj.data(), qQ05per4pie0, f, m2, u);
			UpotXpoles += u;
		}
		// Charge-Dipole
		for (unsigned sj = 0; sj < nd2; sj++) {
			const std::array<double,3> djj = mj.dipole_d_abs(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.dipole_e(sj);
			PotForceChargeDipole(drs, dr2, ejj.data(), minusqmyper4pie0, f, m2, u);
			UpotXpoles += u;
		}
	}
	for (unsigned int si = 0; si < nq1; ++si) {
		const std::array<double,3> dii = mi.quadrupole_d_abs(si);
		const std::array<double,3> eii = mi.quadrupole_e(si);

		// Quadrupole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = mj.charge_d_abs(sj);
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForceChargeQuadrupole(drs, dr2, eii.data(), qQ05per4pie0, f, m1, u);
			UpotXpoles += u;
		}
		// Quadrupole-Quadrupole -------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
			double q2075;
			params >> q2075;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.quadrupole_e(sj);
			PotForce2Quadrupole(drs, dr2, eii.data(), ejj.data(), q2075, f, m1, m2, u);
			UpotXpoles += u;
		}
		// Quadrupole-Dipole -----------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			const std::array<double,3> djj = mj.dipole_d_abs(sj);
			double qmy15;
			params >> qmy15;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.dipole_e(sj);
			//for(unsigned short d=0;d<3;++d) drs[d]=-drs[d]; // avoid that and toggle add/sub below
			PotForceDiQuadrupole(drs, dr2, ejj.data(), eii.data(), qmy15, f, m2, m1, u);
			UpotXpoles += u;
		}
	}
	for (unsigned int si = 0; si < nd1; ++si) {
		const std::array<double,3> dii = mi.dipole_d_abs(si);
		const std::array<double,3> eii = mi.dipole_e(si);
		// Dipole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const std::array<double,3> djj = mj.charge_d_abs(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			PotForceChargeDipole(drs, dr2, eii.data(), minusqmyper4pie0, f, m1, u);
			UpotXpoles += u;
		}
		// Dipole-Quadrupole -----------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			//double drs[3];
			const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
			double myq15;
			params >> myq15;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.quadrupole_e(sj);
			PotForceDiQuadrupole(drs, dr2, eii.data(), ejj.data(), myq15, f, m1, m2, u);
			UpotXpoles += u;
		}
		// Dipole-Dipole ---------------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			//double drs[3];
			const std::array<double,3> djj = mj.dipole_d_abs(sj);
			double my2;
			params >> my2;
			double rffac;
			params >> rffac;
			SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
			const std::array<double,3> ejj = mj.dipole_e(sj);
			PotForce2Dipole(drs, dr2, eii.data(), ejj.data(), my2, rffac, f, m1, m2, u, MyRF);
			UpotXpoles += u;
		}
	}

	// check whether all parameters were used
	mardyn_assert(params.eos());
}

#endif /* POTFORCE_H_ */

