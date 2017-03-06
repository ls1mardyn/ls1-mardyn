#ifndef POTFORCE_H_
#define POTFORCE_H_

/**  @file  */

#include "molecules/Comp2Param.h"
#include "molecules/Molecule.h"
#include "ensemble/PressureGradient.h"
#include "Domain.h"

#include "ensemble/GrandCanonical.h"
#include <iostream>
using namespace std;

//Domain _domain2;

/** @brief Calculates force as spring equivalent for upper layer of upper plate. */
inline void PotForceSpring(const double averageZeroYPos, unsigned long maxID, unsigned long minID, const double springConst, unsigned long molID, double currentYPos, double f[3], unsigned long initStatistics, unsigned long simstep)
{
  if(simstep > 2*initStatistics){
	if (molID <= maxID && molID >= minID){
	    // spring force is of interest for y-direction
	    f[0] = 0.0;
	    f[1] = -springConst*(currentYPos - averageZeroYPos);
	    f[2] = 0.0;
	}
	else{
	    f[0] = 0.0;
	    f[1] = 0.0;
	    f[2] = 0.0;
	}
  }
  else{
	if (simstep > 1.5*initStatistics && molID <= maxID && molID >= minID){
	    // spring force is of interest for y-direction
	    double uniformIncrease = (double)(sqrt((long double)((simstep-1.5*initStatistics)*2.0)/(long double)initStatistics));
	    f[0] = 0.0;
	    f[1] = -springConst*uniformIncrease*(currentYPos - averageZeroYPos);
	    f[2] = 0.0;
	}
	else if (simstep > initStatistics && molID <= maxID && molID >= minID){
	    // spring force is of interest for y-direction
	    f[0] = 0.0;
	    f[1] = -0.5*(currentYPos - averageZeroYPos);
	    f[2] = 0.0;
	}
	else{
	    f[0] = 0.0;
	    f[1] = 0.0;
	    f[2] = 0.0;
	}
  } 
}

/** @brief Calculates force as gravity equivalent. */
inline void PotForceGravity(const int gravitationalDir, const double gravitationalForce, double f[3], unsigned long initStatistics, unsigned long simstep)
{
  if(simstep > 2*initStatistics){
	for(int d = 0; d < 3; d++){
		f[d] = 0.0;
		if(d == gravitationalDir)
			f[d] = gravitationalForce;
	}
  }
  else{
	if (simstep > 1.5*initStatistics){
	    // spring force is of interest for y-direction
	    double uniformIncrease = (double)(sqrt((long double)((simstep-1.5*initStatistics)*2.0)/(long double)initStatistics));
	    for(int d = 0; d < 3; d++){
		f[d] = 0.0;
		if(d == gravitationalDir)
			f[d] = uniformIncrease*gravitationalForce;
	    }
	}
	else if (simstep > initStatistics){
	    for(int d = 0; d < 3; d++){
		f[d] = 0.0;
		if(d == gravitationalDir)
			f[d] = 0.05;
	    }
	}
	else{
	    for(int d = 0; d < 3; d++)
		f[d] = 0.0;
	}
  } 
}

/** @brief Calculates shear force. */
inline void PotForceShear(const double shearYmax, const double shearForce, double currentYPos, double f[3], unsigned long initStatistics, unsigned long simstep, unsigned long rampTime)
{
  double slowAccelerate = 1.0;	
  if((simstep-initStatistics) < rampTime)
			slowAccelerate = (double)(simstep-initStatistics)/rampTime;
  
  //double shearForceTarget = shearForce * shearYmax/2 - fabs((shearYmax/2 - currentYPos) * shearForce);
  //double auxDist = (1-(fabs((shearYmax/2 - currentYPos)/(shearYmax/2))));
  //double shearForceTarget = shearForce * auxDist * auxDist;
//   double shearForceTarget;
//   if (currentYPos < shearYmax/2)
// 	  shearForceTarget = shearForce/10.0 * (-1)/(currentYPos - shearYmax/2 - 0.1);
//   else if (currentYPos > shearYmax/2)
// 	  shearForceTarget = shearForce/10.0 * 1/(currentYPos - shearYmax/2 + 0.1);
//   else
// 	  shearForceTarget = shearForce;
  
  double shearForceTarget = 0.0;
  if (fabs(shearYmax/2 - currentYPos) < 0.5 )
	  shearForceTarget = shearForce;
		  
  shearForceTarget = slowAccelerate*shearForceTarget;
  
//   cout << " shearForce " << shearForceTarget << " y " << round(currentYPos) << endl;
	
  if(simstep > 2*initStatistics){
	    f[0] = shearForceTarget;
	    f[1] = 0.0;
	    f[2] = 0.0;
  }
  else{
	if (simstep > 1.5*initStatistics){
	    // spring force is of interest for y-direction
	    double uniformIncrease = (double)(sqrt((long double)((simstep-1.5*initStatistics)*2.0)/(long double)initStatistics));
	    f[0] = uniformIncrease*shearForceTarget;
	    f[1] = 0.0;
	    f[2] = 0.0;
	}
	else if (simstep > initStatistics){
	    // spring force is of interest for y-direction
	    f[0] = -0.005*shearForceTarget;
	    f[1] = 0.0;
	    f[2] = 0.0;
	}
	else{
	    f[0] = 0.0;
	    f[1] = 0.0;
	    f[2] = 0.0;
	}
  } 
}

/** @brief Calculates force as damper equivalent for upper layer of upper plate.
inline void PotForceDamper(unsigned long maxID, unsigned long minID, const double dampConst, unsigned long molID, double currentVy, double f[3])
{
	if (molID <= maxID && molID >= minID){
	    // spring force is of interest for y-direction
	    f[0] = 0.0;
	    f[1] = -dampConst*currentYy;
	    f[2] = 0.0;
	}
	else{
	    f[0] = 0.0;
	    f[1] = 0.0;
	    f[2] = 0.0;
	}
}	 */

/** @brief Calculates potential and force between 2 Lennard-Jones 12-6 centers. */
inline void PotForceLJ(const double dr[3], const double& invdr2,
                       const double& eps24, const double& sig2,
                       double& sfac, double& ffac, double f[3], double& u6)
{
	double lj6 = sig2 * invdr2; lj6 = lj6 * lj6 * lj6;
	double lj12 = lj6 * lj6;
	double lj12m6 = lj12 - lj6;
	u6 = eps24 * lj12m6;
	
	double fac = eps24 * (lj12 + lj12m6) * invdr2;
	for (unsigned short d = 0; d < 3; ++d)
		f[d] = fac * dr[d];
		
        ffac = eps24 * (lj12 + lj12m6) * invdr2;
	for(unsigned short d = 0; d < 3; d++) f[d] = ffac * dr[d];
        sfac = eps24 * invdr2 * (26.0*lj12 - 7.0*lj6);
}


/** @brief Calculate potential and force between 2 Dipoles. */
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


/** @brief Calculate potential and force between 2 Quadrupoles. */
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


/** @brief Calculate potential and force between a Dipole and Quadrupole. */
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


/** @brief Calculate potential and force between two point charges. */
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


/** @brief Calculate potential and force between an electric charge and a quadrupole. */
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


/** @brief Calculate potential and force between an electric charge and a dipole. */
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
 * @param[in]  calculateLJ    enable or disable calculation of Lennard Jones interactions
 */
inline void PotForce(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, double& Virial, double& VirialIX, double& VirialIY, double& VirialIZ, double& VirialIILL, double& VirialIILM, bool calculateLJ, Domain& domain)
// ???better calc Virial, when molecule forces are calculated:<< endl;
//    summing up molecule virials instead of site virials???
{ // Force Calculation
        double sfac, ffac, tfac, steric, drx2, dry2, drz2;
        
        drx2 = drm[0]*drm[0];
        dry2 = drm[1]*drm[1];
        drz2 = drm[2]*drm[2];
        steric = drx2*dry2 + drx2*drz2 + dry2*drz2;
        
	double f[3];
	double u;
	double drs[3], dr2; // site distance vector & length^2
	long double virialForce;
	long double virialHeat_i, virialHeat_j;
	long double convectivePotHeat_i, convectivePotHeat_j;
	// LJ centers
	// no LJ interaction between solid atoms of the same component
	if ((mi.numTersoff() == 0) || (mi.componentid() != mj.componentid())) {
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
                                        double invdr2 = 1. / dr2;
					PotForceLJ(drs, invdr2, eps24, sig2, sfac, ffac, f, u);  // sfac = (d^2 u / dr^2), ffac = -du / r dr
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
						Upot6LJ += u;
						for (unsigned short d = 0; d < 3; ++d)
						  Virial += drm[d] * f[d];
						VirialIX += drm[0] * f[0];
                                                VirialIY += drm[1] * f[1];
                                                VirialIZ += drm[2] * f[2];
						  
                                                tfac = sfac + ffac;
                                                VirialIILL += drm[0]*drm[0] * (invdr2*drm[0]*drm[0]*tfac - ffac);
                                                VirialIILL += drm[1]*drm[1] * (invdr2*drm[1]*drm[1]*tfac - ffac);
                                                VirialIILL += drm[2]*drm[2] * (invdr2*drm[2]*drm[2]*tfac - ffac);
                                                
                                                VirialIILM += invdr2 * steric * tfac;
						
						if((mi.isHardyStress() || mj.isHardyStress()) && !(domain.getSimstep() % domain.getStressRecordTimeStep())){
						  //calculation of the bond length fraction per unID
						  string stress ("Stress");
						  string weightingFunc = mi.getWeightingFuncStress();
						  mi.calculateHardyIntersection(drm, mj.r(0), mj.r(1), mj.r(2), &domain, stress, weightingFunc);
						  std::map<unsigned long, long double> bondFrac;
						  bondFrac = mi.getBondFractionUNID();
						  for (unsigned short d = 0; d < 3; ++d){
						    for (unsigned short e = 0; e < 3; ++e){
							// the prefactor 0.5 in the virialForce is added because this term is calculated twice (once by molecule mi, once by molecule mj)
							virialForce = 0.5 * drm[d] * f[e];
							virialHeat_i = virialForce;
							virialHeat_j = virialForce;
							for(std::map<unsigned long, long double>::iterator it=bondFrac.begin(); it!=bondFrac.end(); ++it){ 
							    mi.addVirialForceHardyStress(d, e, it->first, virialForce*it->second);
							    mj.addVirialForceHardyStress(d, e, it->first, virialForce*it->second);
							    mi.addDiffusiveHeatfluxHardyStress(d, e, it->first, virialHeat_i*it->second);
							    mj.addDiffusiveHeatfluxHardyStress(d, e, it->first, virialHeat_j*it->second);
							}
						    }
						    convectivePotHeat_i = 0.5 * u/6;
						    convectivePotHeat_j = 0.5 * u/6;
						    mi.addConvectivePotHeatfluxStress(d, convectivePotHeat_i);
						    mj.addConvectivePotHeatfluxStress(d, convectivePotHeat_j);
						  }
						}else if((domain.isStressCalculating(mi.componentid()) || domain.isStressCalculating(mj.componentid())) && !(domain.getSimstep() % domain.getStressRecordTimeStep())){
						  for (unsigned short d = 0; d < 3; ++d){
						    for (unsigned short e = 0; e < 3; ++e){
							virialForce = 0.5 * drm[d] * f[e];
							mi.addVirialForce(d, e, virialForce);
							mj.addVirialForce(d, e, virialForce);
							virialHeat_i = virialForce;// * (mi.v(e) - mi.getDirectedVelocityConfinement(e));
							virialHeat_j = virialForce;// * (mj.v(e) - mj.getDirectedVelocityConfinement(e));
							mi.addDiffusiveHeatfluxStress(d, e, virialHeat_i);
							mj.addDiffusiveHeatfluxStress(d, e, virialHeat_j);
						    }
						    convectivePotHeat_i = 0.5 * u/6;
						    convectivePotHeat_j = 0.5 * u/6;
						    mi.addConvectivePotHeatfluxStress(d, convectivePotHeat_i);
						    mj.addConvectivePotHeatfluxStress(d, convectivePotHeat_j);
						  }
						}
						  
						if(domain.isBulkPressure(mi.componentid()) || domain.isBulkPressure(mj.componentid())){
						  if((mi.r(0) >= domain.getBulkBoundary(0)-domain.getCutoffRadius() && mi.r(0) <= domain.getBulkBoundary(1)+domain.getCutoffRadius() && mi.r(1) >= domain.getBulkBoundary(2)-domain.getCutoffRadius() && mi.r(1) <= domain.getBulkBoundary(3)+domain.getCutoffRadius())
						    || (mj.r(0) >= domain.getBulkBoundary(0)-domain.getCutoffRadius() && mj.r(0) <= domain.getBulkBoundary(1)+domain.getCutoffRadius() && mj.r(1) >= domain.getBulkBoundary(2)-domain.getCutoffRadius() && mj.r(1) <= domain.getBulkBoundary(3)+domain.getCutoffRadius())){
						    for (unsigned short d = 0; d < 3; ++d){
							virialForce = 0.5 * drm[d] * f[d];
							mi.addPressureVirial(d, virialForce);
							mj.addPressureVirial(d, virialForce);
						    }		      
						  }
						}
						
						if((domain.isBarostat(mi.componentid()) || domain.isBarostat(mj.componentid())) && (domain.getSimstep() >= domain.getBarostatTimeInit() && domain.getSimstep() <= domain.getBarostatTimeEnd())){
						  if((mi.r(0) >= domain.getControl_bottom(0) && mi.r(0) <= domain.getControl_top(0) && mi.r(1) >= domain.getControl_bottom(1) && mi.r(1) <= domain.getControl_top(1))
						    || (mj.r(0) >= domain.getControl_bottom(0) && mj.r(0) <= domain.getControl_top(0) && mj.r(1) >= domain.getControl_bottom(1) && mj.r(1) <= domain.getControl_top(1))){
						    for (unsigned short d = 0; d < 3; ++d){
							virialForce = 0.5 * drm[d] * f[d];
							mi.addPressureVirial_barostat(d, virialForce);
							mj.addPressureVirial_barostat(d, virialForce);
						    }
						  }
						}
						
						if((mi.isHardyConfinement() || mj.isHardyConfinement()) && (domain.isConfinement(mi.componentid()) || domain.isConfinement(mj.componentid())) && !(domain.getSimstep() % domain.getConfinementRecordTimeStep())){
						  if((mi.r(1) >= domain.get_confinementMidPoint(3)-domain.getCutoffRadius() && mi.r(1) <= domain.get_confinementMidPoint(1)+domain.getCutoffRadius())
						    || (mj.r(1) >= domain.get_confinementMidPoint(3)-domain.getCutoffRadius() && mj.r(1) <= domain.get_confinementMidPoint(1)+domain.getCutoffRadius())){
						    //calculation of the bond length fraction per unID
						    string stress ("Confinement");
						    string weightingFunc = mi.getWeightingFuncConfinement();
						    if(domain.isConfinementHardy(mi.componentid()) || domain.isConfinementHardy(mj.componentid())){
						      mi.calculateHardyIntersection(drm, mj.r(0), mj.r(1), mj.r(2), &domain, stress, weightingFunc);
						    }
						    std::map<unsigned long, long double> bondFrac;
						    if(domain.isConfinementHardy(mi.componentid()) || domain.isConfinementHardy(mj.componentid()))
						      bondFrac = mi.getBondFractionUNID();
						    for (unsigned short d = 0; d < 3; ++d){
							// the prefactor 0.5 in the virialForce is added because this term is calculated twice (once by molecule mi, once by molecule mj)
							virialForce = 0.5 * drm[d] * f[d];
							mi.addPressureVirialConfinement(d, virialForce);
							mj.addPressureVirialConfinement(d, virialForce);
							if(domain.isConfinementHardy(mi.componentid()) || domain.isConfinementHardy(mj.componentid())){
							  for (unsigned e = 0; e < 3; ++e){
							    virialForce = 0.5 * drm[d] * f[e];
							    virialHeat_i = virialForce;// * (mi.v(e) - mi.getDirectedVelocityConfinement(e));
							    virialHeat_j = virialForce;// * (mj.v(e) - mj.getDirectedVelocityConfinement(e));
							    for(std::map<unsigned long,long double>::iterator it=bondFrac.begin(); it!=bondFrac.end(); ++it){
							      mi.addVirialForceHardyConfinement(d, e, it->first, virialForce*it->second);
							      mj.addVirialForceHardyConfinement(d, e, it->first, virialForce*it->second);
							      mi.addDiffusiveHeatfluxHardyConfinement(d, e, it->first, virialHeat_i*it->second);
							      mj.addDiffusiveHeatfluxHardyConfinement(d, e, it->first, virialHeat_j*it->second);
							    }
							  }
							  convectivePotHeat_i = 0.5 * u/6;// * (mi.v(d) - mi.getDirectedVelocityConfinement(d));
							  convectivePotHeat_j = 0.5 * u/6;// * (mj.v(d) - mj.getDirectedVelocityConfinement(d));
							  mi.addConvectivePotHeatflux(d, convectivePotHeat_i);
							  mj.addConvectivePotHeatflux(d, convectivePotHeat_j);
							}
						    }
						  }		      
						}else if((domain.isConfinement(mi.componentid()) || domain.isConfinement(mj.componentid())) && !(domain.getSimstep() % domain.getConfinementRecordTimeStep())){
						  if((mi.r(1) >= domain.get_confinementMidPoint(3)-domain.getConfinementEdge(5)-domain.getCutoffRadius() && mi.r(1) <= domain.get_confinementMidPoint(1)+domain.getConfinementEdge(5)+domain.getCutoffRadius())
						    || (mj.r(1) >= domain.get_confinementMidPoint(3)-domain.getConfinementEdge(5)-domain.getCutoffRadius() && mj.r(1) <= domain.get_confinementMidPoint(1)+domain.getConfinementEdge(5)+domain.getCutoffRadius())){
						    for (unsigned short d = 0; d < 3; ++d){
							virialForce = 0.5 * drm[d] * f[d];
							mi.addPressureVirialConfinement(d, virialForce);
							mj.addPressureVirialConfinement(d, virialForce);
							for (unsigned e = 0; e < 3; ++e){
							  virialForce = 0.5 * drm[d] * f[e];
							  mi.addVirialForceConfinement(d, e, virialForce);
							  mj.addVirialForceConfinement(d, e, virialForce);
							  virialHeat_i = virialForce;// * (mi.v(e) - mi.getDirectedVelocityConfinement(e));
							  virialHeat_j = virialForce;// * (mj.v(e) - mj.getDirectedVelocityConfinement(e));
							  mi.addDiffusiveHeatflux(d, e, virialHeat_i);
							  mj.addDiffusiveHeatflux(d, e, virialHeat_j);
							}
							convectivePotHeat_i = 0.5 * u/6;// * (mi.v(d) - mi.getDirectedVelocityConfinement(d));
							convectivePotHeat_j = 0.5 * u/6;// * (mj.v(d) - mj.getDirectedVelocityConfinement(d));
							mi.addConvectivePotHeatflux(d, convectivePotHeat_i);
							mj.addConvectivePotHeatflux(d, convectivePotHeat_j);
						    }
						  }		      
						}
					}
				}
			}
		}
	}
	
	//TEST: spring force is added to the LJ-force
	// the molecules to be effected by the spring force lay in one plane and have continous IDs;
	// to add the spring force just once to each molecule of the upper layer of the upper plate
	// an if-request is upstreamed
	if (domain.getPG()->isSpringDamped()){
	    if (mi.getCounter() == 0){
		for (unsigned int si = 0; si < 1; ++si) {
		    PotForceSpring(domain.getPG()->getAverageY(), domain.getPG()->getMaxSpringID(), domain.getPG()->getMinSpringID(), domain.getPG()->getSpringConst(), mi.id(), mi.r(1), f, domain.getInitStatistics(), domain.getSimstep());
		    mi.Fljcenteradd(si, f);
		    for (unsigned short d = 0; d < 3; ++d)
		    {
			Virial += drm[d] * f[d];		// TODO: Check if random or directed virial
			mi.setF_Spring(d, 0.0);
		    }
		
		    for (unsigned short d = 0; d < 3; ++d)
			mi.setF_Spring(d, f[d]);
		    
		    // Prohibits that the spring force is added to one molecule more than once
		    mi.setCounter(1);
		    
		}
	    }
	}
	
	//TEST: gravitational force is added to the LJ-force
	if (domain.getPG()->isGravity()){
	    if (mi.getCounterGravity() == 0 && mi.componentid() == domain.getPG()->getGravityComp()){
		for (unsigned int si = 0; si < 1; ++si) {
		    PotForceGravity(domain.getPG()->getGravityDir(), domain.getPG()->getGravityForce(), f, domain.getInitStatistics(), domain.getSimstep());
		    mi.Fljcenteradd(si, f);
		    for (unsigned short d = 0; d < 3; ++d)
			Virial += drm[d] * f[d];		// TODO: Check if random or directed virial
		    
		    // Prohibits that the gravitational force is added to one molecule more than once
		    mi.setCounterGravity(1);
		    
		}
	    }
	}
	
	//TEST: shear force is added to the LJ-force
	if (domain.getPG()->isShearForce()){
	    if (mi.getCounterShear() == 0){
// 		    cout << " box0 " << domain.getPG()->getShearRateBox(0) << " box1 " << domain.getPG()->getShearRateBox(1) << " box2 " << domain.getPG()->getShearRateBox(2)<< " box3 " << domain.getPG()->getShearRateBox(3) << " width " <<  domain.getPG()->getShearWidth() << endl;
	     if(mi.r(0) > domain.getPG()->getShearRateBox(0) && mi.r(0) < domain.getPG()->getShearRateBox(1) 
		&& mi.r(1) > domain.getPG()->getShearRateBox(2)+domain.getPG()->getShearWidth() && mi.r(1) < domain.getPG()->getShearRateBox(3)-domain.getPG()->getShearWidth()){
		for (unsigned int si = 0; si < 1; ++si) {
// 			cout << " TESTESTEST\n";
		    PotForceShear(domain.getPG()->getShearRateBox(3)-domain.getPG()->getShearRateBox(2)-2*domain.getPG()->getShearWidth(), domain.getPG()->getShearRate(), mi.r(1), f, domain.getInitStatistics(), domain.getSimstep(), domain.getPG()->getShearRampTime());
		    mi.Fljcenteradd(si, f);
		    for (unsigned short d = 0; d < 3; ++d)
			Virial += drm[d] * f[d];		// TODO: Check if random or directed virial
		    
		    // Prohibits that the shear force is added to one molecule more than once
		    mi.setCounterShear(1);
		    
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
                        VirialIX += drm[0] * f[0];
                        VirialIY += drm[1] * f[1];
                        VirialIZ += drm[2] * f[2];
		}
		// Charge-Quadrupole
		for (unsigned sj = 0; sj < nq2; sj++) {
			const double* djj = mj.quadrupole_d(sj);
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.quadrupole_e(sj);
			PotForceChargeQuadrupole(drs, dr2, ejj, qQ05per4pie0, f, m2, u);

			mi.Fchargeadd(si, f);
			mj.Fquadrupolesub(sj, f);
			mj.Madd(m2);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial += drm[d] * f[d];
                        VirialIX += drm[0] * f[0];
                        VirialIY += drm[1] * f[1];
                        VirialIZ += drm[2] * f[2];
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
                        VirialIX += drm[0] * f[0];
                        VirialIY += drm[1] * f[1];
                        VirialIZ += drm[2] * f[2];
		}
	}
	for (unsigned int si = 0; si < nq1; ++si) {
		const double* dii = mi.quadrupole_d(si);
		const double* eii = mi.quadrupole_e(si);

		// Quadrupole-Charge
		for (unsigned sj = 0; sj < ne2; sj++) {
			const double* djj = mj.charge_d(sj);
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			minusSiteSiteDistance(drm, dii, djj, drs, dr2);
			PotForceChargeQuadrupole(drs, dr2, eii, qQ05per4pie0, f, m1, u);

			mi.Fquadrupolesub(si, f);
			mj.Fchargeadd(sj, f);
			mi.Madd(m1);

			UpotXpoles += u;
			for (unsigned short d = 0; d < 3; d++)
				Virial -= drm[d] * f[d];
                        VirialIX -= drm[0] * f[0];
                        VirialIY -= drm[1] * f[1];
                        VirialIZ -= drm[2] * f[2];
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
                        VirialIX += drm[0] * f[0];
                        VirialIY += drm[1] * f[1];
                        VirialIZ += drm[2] * f[2];
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
                        VirialIX -= drm[0] * f[0];
                        VirialIY -= drm[1] * f[1];
                        VirialIZ -= drm[2] * f[2];
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
                        VirialIX -= drm[0] * f[0];
                        VirialIY -= drm[1] * f[1];
                        VirialIZ -= drm[2] * f[2];
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
                        VirialIX += drm[0] * f[0];
                        VirialIY += drm[1] * f[1];
                        VirialIZ += drm[2] * f[2];
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
                        VirialIX += drm[0] * f[0];
                        VirialIY += drm[1] * f[1];
                        VirialIZ += drm[2] * f[2];
		}
	}

	// check whether all parameters were used
	assert(params.eos());
}

inline void PotForce(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, double& Virial, bool calculateLJ, Domain& domain)
{
   double VirIX = 0.0;
   double VirIY = 0.0;
   double VirIZ = 0.0;
   double VirIILL = 0.0;
   double VirIILM = 0.0;
   
   PotForce(mi, mj, params, drm, Upot6LJ, UpotXpoles, MyRF, Virial, VirIX, VirIY, VirIZ, VirIILL, VirIILM, calculateLJ, domain);
   Virial += VirIX + VirIY + VirIZ;
}

/** @brief Calculates the LJ and electrostatic potential energy of the mi-mj interaction (no multi-body potentials are considered) */
inline void FluidPot(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, bool calculateLJ)
{
        double sfac, ffac;
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
				params >> shift6; // must be 0.0 for full LJ

				if (calculateLJ) {
                                        double invdr2 = 1. / dr2;
					PotForceLJ(drs, invdr2, eps24, sig2, sfac, ffac, f, u);
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
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			SiteSiteDistance(drm, dii, djj, drs, dr2);
			const double* ejj = mj.quadrupole_e(sj);
			PotForceChargeQuadrupole(drs, dr2, ejj, qQ05per4pie0, f, m2, u);
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
			double qQ05per4pie0; // 4pie0 = 1 in reduced units
			params >> qQ05per4pie0;
			minusSiteSiteDistance(drm, dii, djj, drs, dr2);
			PotForceChargeQuadrupole(drs, dr2, eii, qQ05per4pie0, f, m1, u);
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

//! @brief Calculate the Tersoff potential for a single atom
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

#endif /* POTFORCE_H_ */

