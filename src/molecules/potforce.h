/***************************************************************************
 *   Copyright (C) 2005 by Martin Bernreuther   *
 *   Martin.Bernreuther@informatik.uni-stuttgart.de   *
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
 ***************************************************************************/

#ifndef POTFORCE_H_
#define POTFORCE_H_

#include "molecules/Comp2Param.h"

/// helper function to calculate the distance between 2 sites
inline void SiteSiteDistance(const double drm[3], const double ds1[3], const double ds2[3], double drs[3], double& dr2)
{
  dr2=0.;
  for(unsigned short d=0;d<3;++d)
  {
    drs[d]=drm[d]+ds1[d]-ds2[d];
    dr2+=drs[d]*drs[d];
  }
}

/// calculate potential and force between 2 Lennard-Jones 12-6 centers
//inline void PotForceLJ(const double dr[3], const double& dr2, ParaStrm& params, double f[3], double& u)
inline void PotForceLJ(const double dr[3], const double& dr2
                      , const double& eps24, const double& sig2
                      ,double f[3], double& u6)
{
  //double eps24;
  //params >> eps24;
  //double sig2;
  //params >> sig2;
  double invdr2=1./dr2;
  double lj6=sig2*invdr2; lj6=lj6*lj6*lj6;
  double lj12=lj6*lj6;
  double lj12m6=lj12-lj6;
  u6=eps24*lj12m6;
  double fac=eps24*(lj12+lj12m6)*invdr2;
  for(unsigned short d=0;d<3;++d) f[d]=fac*dr[d];
}

inline void PotForceLJ(const double dr[3]
                      , const double& eps24, const double& sig2
                      , double f[3], double& u6)
{
  double dr2=0.;
  for(unsigned short d=0;d<3;++d) dr2+=dr[d]*dr[d];
  PotForceLJ(dr,dr2,eps24,sig2,f,u6);
}

/// calculate potential and force between 2 Dipoles (dr2 given)
inline void PotForce2Dipole(const double dr[3], const double& dr2, const double* eii, const double* ejj
                            , const double& my2, const double& rffac
                            , double f[3], double m1[3], double m2[3], double& u, double& MyRF)
{
  double invdr2=1./dr2;
  double invdr1=sqrt(invdr2);
  double myfac=my2*invdr2*invdr1;
  double costi=0.;
  double costj=0.;
  double cosgij=0.;
  for(unsigned short d=0;d<3;++d)
  {
    const double& drd=dr[d];
    const double& eiid=eii[d];
    const double& ejjd=ejj[d];
    costi+=eiid*drd;
    costj+=ejjd*drd;
    cosgij+=eiid*ejjd;
  }
  costi*=invdr1;
  costj*=invdr1;
  u=myfac*(cosgij-3.*costi*costj);
  MyRF-=rffac*cosgij;
  const double partialRijInvdr1=-3.*u*invdr2;
  const double partialTiInvdr1=-myfac*3.*costj*invdr1;
  const double partialTjInvdr1=-myfac*3.*costi*invdr1;
  const double& partialGij=myfac;
  const double fac=-partialRijInvdr1+(costi*partialTiInvdr1+costj*partialTjInvdr1)*invdr1;
  for(unsigned short d=0;d<3;++d) f[d]=fac*dr[d]-partialTiInvdr1*eii[d]-partialTjInvdr1*ejj[d];

  double eiXej[3], eXrij[3];
  eiXej[0]=eii[1]*ejj[2]-eii[2]*ejj[1];
  eiXej[1]=eii[2]*ejj[0]-eii[0]*ejj[2];
  eiXej[2]=eii[0]*ejj[1]-eii[1]*ejj[0];
  eXrij[0]=eii[1]*dr[2]-eii[2]*dr[1];
  eXrij[1]=eii[2]*dr[0]-eii[0]*dr[2];
  eXrij[2]=eii[0]*dr[1]-eii[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m1[d]=-partialTiInvdr1*eXrij[d]+(-partialGij+rffac)*eiXej[d];
  eXrij[0]=ejj[1]*dr[2]-ejj[2]*dr[1];
  eXrij[1]=ejj[2]*dr[0]-ejj[0]*dr[2];
  eXrij[2]=ejj[0]*dr[1]-ejj[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m2[d]=-partialTjInvdr1*eXrij[d]+(partialGij-rffac)*eiXej[d];
}

/// calculate potential and force between 2 Dipoles
inline void PotForce2Dipole(const double dr[3], const double* eii, const double* ejj
                            , const double& my2, const double& rffac
                            ,double f[3], double m1[3], double m2[3], double& u, double& MyRF)
{
  double dr2=0.;
  for(unsigned short d=0;d<3;++d) dr2+=dr[d]*dr[d];
  PotForce2Dipole(dr,dr2,eii,ejj,my2,rffac,f,m1,m2,u,MyRF);
}



/// calculate potential and force between 2 Quadrupoles (dr2 given)
inline void PotForce2Quadrupole(const double dr[3], const double& dr2, const double* eii, const double* ejj
                               , const double& q2075
                               , double f[3], double m1[3], double m2[3], double& u)
{
  double invdr2=1./dr2;
  double invdr1=sqrt(invdr2);
  double qfac=q2075*invdr2*invdr2*invdr1;
  double costi=0.;
  double costj=0.;
  double cosgij=0.;
  for(unsigned short d=0;d<3;++d)
  {
    const double& drd=dr[d];
    const double& eiid=eii[d];
    const double& ejjd=ejj[d];
    costi+=eiid*drd;
    costj+=ejjd*drd;
    cosgij+=eiid*ejjd;
  }
  costi*=invdr1;
  costj*=invdr1;
  double cos2ti=costi*costi;
  double cos2tj=costj*costj;
  double term=(cosgij-5.*costi*costj);
  u=qfac*(1.-5.*(cos2ti+cos2tj)-15.*cos2ti*cos2tj+2.*term*term);
  const double partialRijInvdr1=-5.*u*invdr2;
  const double partialTiInvdr1=-qfac*10.*(costi+3.*costi*cos2tj+2.*costj*term)*invdr1;
  const double partialTjInvdr1=-qfac*10.*(costj+3.*cos2ti*costj+2.*costi*term)*invdr1;
  const double& partialGij=qfac*4.*term;
  const double fac=-partialRijInvdr1+(costi*partialTiInvdr1+costj*partialTjInvdr1)*invdr1;
  for(unsigned short d=0;d<3;++d) f[d]=fac*dr[d]-partialTiInvdr1*eii[d]-partialTjInvdr1*ejj[d];

  double eiXej[3], eXrij[3];
  eiXej[0]=eii[1]*ejj[2]-eii[2]*ejj[1];
  eiXej[1]=eii[2]*ejj[0]-eii[0]*ejj[2];
  eiXej[2]=eii[0]*ejj[1]-eii[1]*ejj[0];
  eXrij[0]=eii[1]*dr[2]-eii[2]*dr[1];
  eXrij[1]=eii[2]*dr[0]-eii[0]*dr[2];
  eXrij[2]=eii[0]*dr[1]-eii[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m1[d]=-partialTiInvdr1*eXrij[d]-partialGij*eiXej[d];

  eXrij[0]=ejj[1]*dr[2]-ejj[2]*dr[1];
  eXrij[1]=ejj[2]*dr[0]-ejj[0]*dr[2];
  eXrij[2]=ejj[0]*dr[1]-ejj[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m2[d]=-partialTjInvdr1*eXrij[d]+partialGij*eiXej[d];
}

/// calculate potential and force between 2 Quadrupoles
inline void PotForce2Quadrupole(const double dr[3], const double* eii, const double* ejj
                               , const double& q2075
                               , double f[3], double m1[3], double m2[3], double& u)
{
  double dr2=0.;
  for(unsigned short d=0;d<3;++d) dr2+=dr[d]*dr[d];
  PotForce2Quadrupole(dr,dr2,eii,ejj,q2075,f,m1,m2,u);
}


/// calculate potential and force between a Dipole and Quadrupole (dr2 given)
inline void PotForceDiQuadrupole(const double dr[3], const double& dr2, const double* eii, const double* ejj
                                , const double& myq15
                                , double f[3], double m1[3], double m2[3], double& u)
{
  double invdr2=1./dr2;
  double invdr1=sqrt(invdr2);
  double myqfac=myq15*invdr2*invdr2;
  double costi=0.;
  double costj=0.;
  double cosgij=0.;
  for(unsigned short d=0;d<3;++d)
  {
    const double& drd=dr[d];
    const double& eiid=eii[d];
    const double& ejjd=ejj[d];
    costi+=eiid*drd;
    costj+=ejjd*drd;
    cosgij+=eiid*ejjd;
  }
  costi*=invdr1;
  costj*=invdr1;
  double cos2tj=costj*costj;
  u=myqfac*(-costi*(5.*cos2tj-1.)+2.*cosgij*costj);
  const double partialRijInvdr1=-4.*u*invdr2;
  const double partialTiInvdr1=myqfac*(-5.*cos2tj+1.)*invdr1;
  const double partialTjInvdr1=myqfac*2.*(-5.*costi*costj+cosgij)*invdr1;
  const double& partialGij=myqfac*2.*costj;
  const double fac=-partialRijInvdr1+(costi*partialTiInvdr1+costj*partialTjInvdr1)*invdr1;
  for(unsigned short d=0;d<3;++d) f[d]=fac*dr[d]-partialTiInvdr1*eii[d]-partialTjInvdr1*ejj[d];

  double eiXej[3], eXrij[3];
  eiXej[0]=eii[1]*ejj[2]-eii[2]*ejj[1];
  eiXej[1]=eii[2]*ejj[0]-eii[0]*ejj[2];
  eiXej[2]=eii[0]*ejj[1]-eii[1]*ejj[0];
  eXrij[0]=eii[1]*dr[2]-eii[2]*dr[1];
  eXrij[1]=eii[2]*dr[0]-eii[0]*dr[2];
  eXrij[2]=eii[0]*dr[1]-eii[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m1[d]=-partialTiInvdr1*eXrij[d]-partialGij*eiXej[d];

  eXrij[0]=ejj[1]*dr[2]-ejj[2]*dr[1];
  eXrij[1]=ejj[2]*dr[0]-ejj[0]*dr[2];
  eXrij[2]=ejj[0]*dr[1]-ejj[1]*dr[0];
  for(unsigned short d=0;d<3;++d) m2[d]=-partialTjInvdr1*eXrij[d]+partialGij*eiXej[d];
}

/// calculate potential and force between a Dipole and Quadrupole
inline void PotForceDiQuadrupole(const double dr[3], const double* eii, const double* ejj
                                 , const double& myq15
                                 , double f[3], double m1[3], double m2[3], double& u)
{
  double dr2=0.;
  for(unsigned short d=0;d<3;++d) dr2+=dr[d]*dr[d];
  PotForce2Quadrupole(dr,dr2,eii,ejj,myq15,f,m1,m2,u);
}



/** calculate Potential and Force between molecules (all site-site interactions)
    paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
    (cmp. comp2param.h/.cpp)
*/
inline void PotForce(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, double& Virial )
// ???better calc Virial, when molecule forces are calculated:
//    summing up molecule virials instead of site virials???
{ // Force Calculation
  double f[3];
  double u;
  double drs[3],dr2;  // site distance vector & length^2
  // LJ centers
  const unsigned int nc1=mi.numLJcenters();
  const unsigned int nc2=mj.numLJcenters();
  for(unsigned int si=0;si<nc1;++si)
  {
    const double* dii=mi.ljcenter_d(si);
    for(unsigned int sj=0;sj<nc2;++sj)
    {
      const double* djj=mj.ljcenter_d(sj);
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      double eps24;
      params >> eps24;
      double sig2;
      params >> sig2;
      PotForceLJ(drs,dr2,eps24,sig2,f,u);
// even for interactions within the cell a neighbor might try to add/subtract
// better use atomic...
// and even better use a order where critical sections occure only at some boundary cells...
#ifdef _OPENMP
#pragma omp critical
#endif
{
      mi.Fljcenteradd(si,f);
      mj.Fljcentersub(sj,f);}
      Upot6LJ+=u;
      /*
      u/=6.;
      mi.Upotadd(u);
      mj.Upotadd(u);
      */
      for(unsigned short d=0;d<3;++d) Virial+=drm[d]*f[d];
    }
  }

  double m1[3],m2[3]; // moments

  const unsigned int nd1=mi.numDipoles();
  const unsigned int nd2=mj.numDipoles();
  const unsigned int nq1=mi.numQuadrupoles();
  const unsigned int nq2=mj.numQuadrupoles();
  for(unsigned int si=0;si<nd1;++si)
  {
    const double* dii=mi.dipole_d(si);
    // Dipole-Dipole ---------------------------
    for(unsigned int sj=0;sj<nd2;++sj)
    {
      //double drs[3];
      const double* djj=mj.dipole_d(sj);
      double my2;
      params >> my2;
      double rffac;
      params >> rffac;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* eii=mi.dipole_e(si);
      const double* ejj=mj.dipole_e(sj);
      PotForce2Dipole(drs,dr2,eii,ejj,my2,rffac,f,m1,m2,u,MyRF);

      mi.Fdipoleadd(si,f);
      mj.Fdipolesub(sj,f);
      mi.Madd(m1);
      mj.Madd(m2);
      UpotXpoles+=u;
      /*
      mi.Upotadd(u);
      mj.Upotadd(u);
      */
      for(unsigned short d=0;d<3;++d) Virial+=drm[d]*f[d];
    }
    // Dipole-Quadrupole -----------------------
    for(unsigned int sj=0;sj<nq2;++sj)
    {
      //double drs[3];
      const double* djj=mj.quadrupole_d(sj);
      double myq15;
      params >> myq15;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* eii=mi.dipole_e(si);
      const double* ejj=mj.quadrupole_e(sj);
      PotForceDiQuadrupole(drs,dr2,eii,ejj,myq15,f,m1,m2,u);

      mi.Fdipoleadd(si,f);
      mj.Fquadrupolesub(sj,f);
      mi.Madd(m1);
      mj.Madd(m2);
      UpotXpoles+=u;
      /*
      mi.Upotadd(u);
      mj.Upotadd(u);
      */
      for(unsigned short d=0;d<3;++d) Virial+=drm[d]*f[d];
    }
  }
  for(unsigned int si=0;si<nq1;++si)
  {
    const double* dii=mi.quadrupole_d(si);
    // Quadrupole-Dipole -----------------------
    for(unsigned int sj=0;sj<nd2;++sj)
    {
      //double drs[3];
      const double* djj=mj.dipole_d(sj);
      double qmy15;
      params >> qmy15;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* eii=mi.quadrupole_e(si);
      const double* ejj=mj.dipole_e(sj);
      //for(unsigned short d=0;d<3;++d) drs[d]=-drs[d]; // avoid that and toggle add/sub below
      PotForceDiQuadrupole(drs,dr2,ejj,eii,qmy15,f,m2,m1,u);

      mi.Fquadrupolesub(si,f);
      mj.Fdipoleadd(sj,f);
      mi.Msub(m1);
      mj.Msub(m2);
      UpotXpoles-=u;
      /*
      mi.Upotsub(u);
      mj.Upotsub(u);
      */
      for(unsigned short d=0;d<3;++d) Virial-=drm[d]*f[d];
    }
    // Quadrupole-Quadrupole -------------------
    for(unsigned int sj=0;sj<nq2;++sj)
    {
      //double drs[3];
      const double* djj=mj.quadrupole_d(sj);
      double q2075;
      params >> q2075;
      SiteSiteDistance(drm,dii,djj,drs,dr2);
      const double* eii=mi.quadrupole_e(si);
      const double* ejj=mj.quadrupole_e(sj);
      PotForce2Quadrupole(drs,dr2,eii,ejj,q2075,f,m1,m2,u);

      mi.Fquadrupoleadd(si,f);
      mj.Fquadrupolesub(sj,f);
      mi.Madd(m1);
      mj.Madd(m2);

      UpotXpoles+=u;
      /*
      mi.Upotadd(u);
      mj.Upotadd(u);
      */
      for(unsigned short d=0;d<3;++d) Virial+=drm[d]*f[d];
    }
  }
  // check, if all parameters were used
  assert(params.eos());
}

//XXX benötigt?
//inline void PotForce(Molecule& mi, Molecule& mj, Comp2Param& comp2params, double L[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, double& Virial)
//{
//  ParaStrm& params=comp2params(mi.componentid(),mj.componentid());
//  params.reset_read();
//  double drm[3];
//  mj.dist2(mi,L,drm);
//  PotForce(mi, mj, params,drm,Upot6LJ,UpotXpoles,MyRF,Virial);
//}

#endif /*POTFORCE_H_*/
