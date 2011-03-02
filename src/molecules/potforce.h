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

#ifndef POTFORCE_H_
#define POTFORCE_H_

#include "molecules/Molecule.h"
#include "molecules/Comp2Param.h"

/// helper function to calculate the distance between 2 sites
inline void SiteSiteDistance(const double drm[3], const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
	for (unsigned short d = 0; d < 3; ++d)
		drs[d] = drm[d] + ds1[d] - ds2[d];
	dr2 = drs[0]*drs[0] + drs[1]*drs[1] + drs[2]*drs[2];
}
inline void minusSiteSiteDistance(const double drm[3], const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
	for (unsigned short d = 0; d < 3; ++d)
		drs[d] = ds2[d] - drm[d] - ds1[d];
	dr2 = drs[0]*drs[0] + drs[1]*drs[1] + drs[2]*drs[2];
}

/// calculate potential and force between 2 Lennard-Jones 12-6 centers
//inline void PotForceLJ(const double dr[3], const double& dr2, ParaStrm& params, double f[3], double& u)
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

inline void PotForceLJ(const double dr[3],
                       const double& eps24, const double& sig2,
                       double f[3], double& u6)
{
	double dr2 = 0.;
	for (unsigned short d = 0; d < 3; ++d)
		dr2 += dr[d] * dr[d];
	PotForceLJ(dr, dr2, eps24, sig2, f, u6);
}

/// calculate potential and force between 2 Dipoles (dr2 given)
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

/// calculate potential and force between 2 Dipoles
inline void PotForce2Dipole(const double dr[3], const double* eii, const double* ejj,
                            const double& my2, const double& rffac,
                            double f[3], double m1[3], double m2[3], double& u, double& MyRF)
{
	double dr2 = 0.;
	for (unsigned short d = 0; d < 3; ++d)
		dr2 += dr[d] * dr[d];
	PotForce2Dipole(dr, dr2, eii, ejj, my2, rffac, f, m1, m2, u, MyRF);
}

/// calculate potential and force between 2 Quadrupoles (dr2 given)
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
	const double& partialGij = qfac * 4. * term;
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

/// calculate potential and force between 2 Quadrupoles
inline void PotForce2Quadrupole(const double dr[3], const double* eii, const double* ejj,
                                const double& q2075,
                                double f[3], double m1[3], double m2[3], double& u)
{
	double dr2 = 0.;
	for (unsigned short d = 0; d < 3; ++d)
		dr2 += dr[d] * dr[d];
	PotForce2Quadrupole(dr, dr2, eii, ejj, q2075, f, m1, m2, u);
}

/// calculate potential and force between a Dipole and Quadrupole (dr2 given)
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
	const double& partialGij = myqfac * 2. * costj;
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

/// calculate potential and force between a Dipole and Quadrupole
inline void PotForceDiQuadrupole(const double dr[3], const double* eii, const double* ejj,
                                 const double& myq15,
                                 double f[3], double m1[3], double m2[3], double& u)
{
	double dr2 = 0.;
	for (unsigned short d = 0; d < 3; ++d)
		dr2 += dr[d] * dr[d];
	PotForce2Quadrupole(dr, dr2, eii, ejj, myq15, f, m1, m2, u);
}

/// calculate potential and force between two point charges (dr2 given)
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

/// calculate potential and force between an electric charge and a quadrupole
inline void PotForceChargeQuadrupole(const double dr[3], const double& dr2,
                                     const double* ejj, const double& qQ025per4pie0,
                                     double f[3], double m2[3], double& u)
{
	double invdr2 = 1.0 / dr2;
	double invdr = sqrt(invdr2);
	double costj = 0;
	for (unsigned short d = 0; d < 3; d++) {
		costj += ejj[d] * dr[d];
	}
	costj *= invdr;
	double qQinv4dr3 = qQ025per4pie0 * invdr * invdr2;
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

/// calculate potential and force between an electric charge and a dipole
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
	const double fac = 3.0 * u * invdr*invdr2;
	for (unsigned short d = 0; d < 3; ++d)
		f[d] = fac*dr[d] - partialTjInvdr1 * ejj[d];

	double minuseXrij[3];
	minuseXrij[0] = ejj[2]*dr[1] - ejj[1]*dr[2];
	minuseXrij[1] = ejj[0]*dr[2] - ejj[2]*dr[0];
	minuseXrij[2] = ejj[1]*dr[0] - ejj[0]*dr[1];
	for (unsigned short d = 0; d < 3; ++d)
		m2[d] = partialTjInvdr1 * minuseXrij[d];
}

/** calculate Potential and Force between molecules (all site-site interactions)
    paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
    (cmp. comp2param.h/.cpp)

    drm == distance FROM j TO i ... !!!
*/
inline void PotForce(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, double& Virial, bool calculateLJ)
// ???better calc Virial, when molecule forces are calculated:
//    summing up molecule virials instead of site virials???
{ // Force Calculation
	double f[3];
	double u;
	double drs[3], dr2; // site distance vector & length^2
	// LJ centers
	// no LJ interaction between solid atoms of the same component
	const unsigned int nt1 = mi.numTersoff();
	if ((mi.componentid() != mj.componentid()) || !nt1) {
		const unsigned int nc1 = mi.numLJcenters();
		const unsigned int nc2 = mj.numLJcenters();
		for (unsigned int si = 0; si < nc1; ++si) {
			const double* dii = mi.ljcenter_d(si);
			for (unsigned int sj = 0; sj < nc2; ++sj) {
				const double* djj = mj.ljcenter_d(sj);
				SiteSiteDistance(drm, dii, djj, drs, dr2);
				double eps24;
				params >> eps24;
				double sig2;
				params >> sig2;
				double shift6;
				params >> shift6; // must be 0.0 for full LJ
				if (calculateLJ) {
					PotForceLJ(drs, dr2, eps24, sig2, f, u);
					u += shift6;

// even for interactions within the cell a neighbor might try to add/subtract
// better use atomic...
// and even better use a order where critical sections occure only at some boundary cells...
#if defined(ENABLE_OPENMP)
#pragma omp critical
#endif
					{
						mi.Fljcenteradd(si, f);
						mj.Fljcentersub(sj, f);
					}
					Upot6LJ += u;
					for (unsigned short d = 0; d < 3; ++d)
						Virial += drm[d] * f[d];
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
		const double* dii = mi.charge_d(si);
		// Charge-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const double* djj = mj.charge_d(sj);
			double q1q2per4pie0; // 4pie0 = 1 in reduced units
			params >> q1q2per4pie0;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);

			mi.Fchargeadd(si, f);
			mj.Fchargesub(sj, f);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial += drm[d] * f[d];
		}
		// Charge-Quadrupole
		for (unsigned sj = 0; sj < nq2; sj++) {
			const double* djj = mj.quadrupole_d(sj);
			double qQ025per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ025per4pie0;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.quadrupole_e(sj);
			PotForceChargeQuadrupole(drs, dr2, ejj, qQ025per4pie0, f, m2, u);

			mi.Fchargeadd(si, f);
			mj.Fquadrupolesub(sj, f);
			mj.Madd(m2);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial += drm[d] * f[d];
		}
		// Charge-Dipole
		for (unsigned sj = 0; sj < nd2; sj++) {
			const double* djj = mj.dipole_d(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.dipole_e(sj);
			PotForceChargeDipole(drs, dr2, ejj, minusqmyper4pie0, f, m2, u);

			mi.Fchargeadd(si, f);
			mj.Fdipolesub(sj, f);
			mj.Madd(m2);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial += drm[d] * f[d];
		}
	}
	for (unsigned int si = 0; si < nq1; ++si) {
		const double* dii = mi.quadrupole_d(si);
		const double* eii = mi.quadrupole_e(si);

		// Quadrupole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const double* djj = mj.charge_d(sj);
			double qQ025per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ025per4pie0;
			minusSiteSiteDistance(drm, dii, djj, drs, dr2);
			PotForceChargeQuadrupole(drs, dr2, eii, qQ025per4pie0, f, m1, u);

			mi.Fquadrupolesub(si, f);
			mj.Fchargeadd(sj, f);
			mi.Madd(m1);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial -= drm[d] * f[d];
		}
		// Quadrupole-Quadrupole -------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			//double drs[3];
			const double* djj = mj.quadrupole_d(sj);
			double q2075;
			params >> q2075;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.quadrupole_e(sj);
			PotForce2Quadrupole(drs, dr2, eii, ejj, q2075, f, m1, m2, u);

			mi.Fquadrupoleadd(si, f);
			mj.Fquadrupolesub(sj, f);
			mi.Madd(m1);
			mj.Madd(m2);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				Virial += drm[d] * f[d];
		}
		// Quadrupole-Dipole -----------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			//double drs[3];
			const double* djj = mj.dipole_d(sj);
			double qmy15;
			params >> qmy15;
			minusSiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.dipole_e(sj);
			PotForceDiQuadrupole(drs, dr2, ejj, eii, qmy15, f, m2, m1, u);

			mi.Fquadrupolesub(si, f);
			mj.Fdipoleadd(sj, f);
			mi.Madd(m1);
			mj.Madd(m2);
			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial -= drm[d] * f[d];
		}
	}
	for (unsigned int si = 0; si < nd1; ++si) {
		const double* dii = mi.dipole_d(si);
		const double* eii = mi.dipole_e(si);
		// Dipole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const double* djj = mj.charge_d(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			minusSiteSiteDistance(drm, dii, djj, drs, dr2);
			PotForceChargeDipole(drs, dr2, eii, minusqmyper4pie0, f, m1, u);

			mi.Fdipolesub(si, f);
			mj.Fchargeadd(sj, f);
			mi.Madd(m1);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial -= drm[d] * f[d];
		}
		// Dipole-Quadrupole -----------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			//double drs[3];
			const double* djj = mj.quadrupole_d(sj);
			double myq15;
			params >> myq15;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.quadrupole_e(sj);
			PotForceDiQuadrupole(drs, dr2, eii, ejj, myq15, f, m1, m2, u);

			mi.Fdipoleadd(si, f);
			mj.Fquadrupolesub(sj, f);
			mi.Madd(m1);
			mj.Madd(m2);
			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				Virial += drm[d] * f[d];
		}
		// Dipole-Dipole ---------------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			const double* djj = mj.dipole_d(sj);
			double my2;
			params >> my2;
			double rffac;
			params >> rffac;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.dipole_e(sj);
			PotForce2Dipole(drs, dr2, eii, ejj, my2, rffac, f, m1, m2, u, MyRF);

			mi.Fdipoleadd(si, f);
			mj.Fdipolesub(sj, f);
			mi.Madd(m1);
			mj.Madd(m2);
			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; ++d)
				Virial += drm[d] * f[d];
		}
	}

	// check whether all parameters were used
	assert(params.eos());
}

/*
 * calculates the LJ and electrostatic potential energy of the mi-mj interaction (no multi-body potentials are considered)
 */
inline void FluidPot(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, bool calculateLJ)
{
	double f[3];
	double u;
	double drs[3], dr2; // site distance vector & length^2
	// no LJ interaction between equal solid atoms
	const unsigned int nt1 = mi.numTersoff();
	if ((mi.componentid() != mj.componentid()) || !nt1) {
		const unsigned int nc1 = mi.numLJcenters();
		const unsigned int nc2 = mj.numLJcenters();
		for (unsigned int si = 0; si < nc1; ++si) {
			const double* dii = mi.ljcenter_d(si);
			for (unsigned int sj = 0; sj < nc2; ++sj) {
				const double* djj = mj.ljcenter_d(sj);
				SiteSiteDistance(drm, dii, djj, drs, dr2);
				double eps24;
				params >> eps24;
				double sig2;
				params >> sig2;
				double shift6;
				/* FIXME: alternative? */
				params >> shift6; // must be 0.0 for full LJ

				if (calculateLJ) {
					PotForceLJ(drs, dr2, eps24, sig2, f, u);
					u += shift6;
					Upot6LJ += u;
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
		const double* dii = mi.charge_d(si);
		// Charge-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const double* djj = mj.charge_d(sj);
			double q1q2per4pie0; // 4pie0 = 1 in reduced units
			params >> q1q2per4pie0;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);
			UpotXpoles += u;
		}
		// Charge-Quadrupole
		for (unsigned sj = 0; sj < nq2; sj++) {
			const double* djj = mj.quadrupole_d(sj);
			double qQ025per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ025per4pie0;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.quadrupole_e(sj);
			PotForceChargeQuadrupole(drs, dr2, ejj, qQ025per4pie0, f, m2, u);
			UpotXpoles += u;
		}
		// Charge-Dipole
		for (unsigned sj = 0; sj < nd2; sj++) {
			const double* djj = mj.dipole_d(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.dipole_e(sj);
			PotForceChargeDipole(drs, dr2, ejj, minusqmyper4pie0, f, m2, u);
			UpotXpoles += u;
		}
	}
	for (unsigned int si = 0; si < nq1; ++si) {
		const double* dii = mi.quadrupole_d(si);
		const double* eii = mi.quadrupole_e(si);

		// Quadrupole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const double* djj = mj.charge_d(sj);
			double qQ025per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ025per4pie0;
			minusSiteSiteDistance(drm, dii, djj, drs, dr2);
			PotForceChargeQuadrupole(drs, dr2, eii, qQ025per4pie0, f, m1, u);
			UpotXpoles += u;
		}
		// Quadrupole-Quadrupole -------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			const double* djj = mj.quadrupole_d(sj);
			double q2075;
			params >> q2075;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.quadrupole_e(sj);
			PotForce2Quadrupole(drs, dr2, eii, ejj, q2075, f, m1, m2, u);
			UpotXpoles += u;
		}
		// Quadrupole-Dipole -----------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			const double* djj = mj.dipole_d(sj);
			double qmy15;
			params >> qmy15;
			minusSiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.dipole_e(sj);
			//for(unsigned short d=0;d<3;++d) drs[d]=-drs[d]; // avoid that and toggle add/sub below
			PotForceDiQuadrupole(drs, dr2, ejj, eii, qmy15, f, m2, m1, u);
			UpotXpoles += u;
		}
	}
	for (unsigned int si = 0; si < nd1; ++si) {
		const double* dii = mi.dipole_d(si);
		const double* eii = mi.dipole_e(si);
		// Dipole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const double* djj = mj.charge_d(sj);
			double minusqmyper4pie0;
			params >> minusqmyper4pie0;
			minusSiteSiteDistance(drm, dii, djj, drs, dr2);
			PotForceChargeDipole(drs, dr2, eii, minusqmyper4pie0, f, m1, u);
			UpotXpoles += u;
		}
		// Dipole-Quadrupole -----------------------
		for (unsigned int sj = 0; sj < nq2; ++sj) {
			//double drs[3];
			const double* djj = mj.quadrupole_d(sj);
			double myq15;
			params >> myq15;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.quadrupole_e(sj);
			PotForceDiQuadrupole(drs, dr2, eii, ejj, myq15, f, m1, m2, u);
			UpotXpoles += u;
		}
		// Dipole-Dipole ---------------------------
		for (unsigned int sj = 0; sj < nd2; ++sj) {
			//double drs[3];
			const double* djj = mj.dipole_d(sj);
			double my2;
			params >> my2;
			double rffac;
			params >> rffac;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.dipole_e(sj);
			PotForce2Dipole(drs, dr2, eii, ejj, my2, rffac, f, m1, m2, u, MyRF);
			UpotXpoles += u;
		}
	}

	// check whether all parameters were used
	assert(params.eos());
}

//!
//! drij should contain the distance from j to i.
//!
inline double Tersoffbij(Molecule* mi, Molecule* mj,
                         double params[15], double drij[3], double drij1)
{
	double drik[3];
	double zeta = 0.0;
	unsigned i_curTN = mi->getCurTN();
	Molecule* mk;

	for (unsigned nmk = 0; nmk < i_curTN; nmk++) {
		mk = mi->getTersoffNeighbour(nmk);
		if ((mk->id() == mj->id()) || (mk->numTersoff() == 0)) continue;

		double drik2 = mk->dist2(*mi, drik); // distance k->i
		if (params[12] > drik2) {
			double drik1 = sqrt(drik2);
			double drdr = 0.0;
			for (int d = 0; d < 3; d++)
				drdr += drij[d] * drik[d];
			double h_min_costheta = params[2] - drdr/(drij1*drik1);
			double g = params[11] - params[3] / (params[4] + h_min_costheta*h_min_costheta);
			if (drik1 > params[0]) {
				double kC = 0.5 + 0.5 * cos((drik1 - params[0]) * params[10]);
				zeta += kC * g;
			}
			else
				zeta += g;
		}
	}
	double BZ = params[8] * zeta;
	double BZtoN = pow(BZ, params[9]);
	double b = pow(1.0 + BZtoN, params[14]);
	return b;
}

//! @brief returns the sum, whereas Uij is added to the double reference
//!
//! drij should contain the distance from j to i.
//!
inline double TersoffUIJplusUJI(Molecule* mi, Molecule* mj, double params[15],
                                double drij[3], double drij2, double& UpotTersoff)
{
	double Uij = 0.0;
	double Uji = 0.0;
	double drji[3];
	double drij1 = sqrt(drij2);
	for (int d = 0; d < 3; d++)
		drji[d] = -drij[d]; // drji: distance i->j

	double UA = params[13] * exp(params[7] * drij1);
	double UR = params[5] * exp(params[6] * drij1);
	double bij = Tersoffbij(mi, mj, params, drij, drij1);
	double bji = Tersoffbij(mj, mi, params, drji, drij1);
	if (drij1 > params[0]) {
		double kChalf = 0.25 + 0.25*cos((drij1 - params[0]) * params[10]);
		Uij = kChalf * (UR + bij*UA);
		Uji = kChalf * (UR + bji*UA);
	}
	else {
		Uij = 0.5 * (UR + bij*UA);
		Uji = 0.5 * (UR + bji*UA);
	}
	UpotTersoff += Uij;
	return Uij + Uji;
}
//!
//! drij should contain the distance from j to i.
//!
inline double TersoffUIJattr(Molecule* mi, Molecule* mj,
                             double params[15], double drij[3], double drij2)
{
	double drij1 = sqrt(drij2);
	double UA = params[13] * exp(params[7] * drij1);
	double bij = Tersoffbij(mi, mj, params, drij, drij1);

	if (drij1 > params[0]) {
		double kChalf = 0.25 + 0.25 * cos((drij1 - params[0]) * params[10]);
		return kChalf * bij * UA;
	}
	else
		return 0.5 * bij * UA;
}

//! @brief calculate Tersoff potential for a single atom
//!
//! parameters: R, S, h, c^2, d^2, A, -lambda, -mu, beta, n_i, pi/(S-R), 1+(c/d)^2, S^2, -B, -0.5/n_i
//! A "Molecule" may have at most a single Tersoff centre.
//!
//! the function itself returns the sum of all Uij, Uji, and the attractive
//! part of Ujk, i.e. all potential energy terms that are influenced by atom i. 
//! However, only the sum over the Uij is added to UpotTersoff.
//!
//! used for computing the Tersoff potential based force on the atom
//!
inline double TersoffPotential(Molecule* mi, double params[15], double& UpotTersoff)
{
	double Ui = 0.0;

	double distanceVector[3];
	unsigned i_curTN = mi->getCurTN();
	Molecule *mj, *mk;

	for (unsigned nmj = 0; nmj < i_curTN; nmj++) {
		mj = mi->getTersoffNeighbour(nmj);
		if (mj->numTersoff() == 0)
			continue;
		double drij2 = mj->dist2(*mi, distanceVector); // distance j->i
		if (params[12] > drij2)
			Ui += TersoffUIJplusUJI(mi, mj, params, distanceVector, drij2, UpotTersoff);

		unsigned j_curTN = mj->getCurTN();

		for (unsigned nmk = 0; nmk < j_curTN; nmk++) {
			mk = mj->getTersoffNeighbour(nmk);
			if (mk->id() == mi->id() || (mk->numTersoff() == 0))
				continue;
			double drjk2 = mk->dist2(*mj, distanceVector); // distance k->j
			if (params[12] > drjk2)
				Ui += TersoffUIJattr(mj, mk, params, distanceVector, drjk2);
		}
	}

	return Ui;
}
inline void TersoffPotForce(Molecule* mi, double params[15], double& UpotTersoff, double delta_r)
{
	double f[3];
	if (mi->numTersoff() == 0) return;
	double offsetU = TersoffPotential(mi, params, UpotTersoff);

	for (int d = 0; d < 3; d++) {
		mi->move(d, delta_r);

		double irrelevantU = 0.0;
		double currentU = TersoffPotential(mi, params, irrelevantU);

		f[d] = (offsetU - currentU) / delta_r;
		mi->move(d, -delta_r);
	}

	mi->Ftersoffadd(0, f);
}

#endif /*POTFORCE_H_*/

