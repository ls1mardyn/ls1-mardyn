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
class Molecule{

public:
  enum streamtype { RESTART };

  Molecule(unsigned long id=0, int componentid=0
          , double rx=0., double ry=0., double rz=0.
          , double vx=0., double vy=0., double vz=0.
          , double q0=0., double q1=0., double q2=0., double q3=0.
          , double Dx=0., double Dy=0., double Dz=0.
          , const std::vector<Component>* components=NULL
          );
  Molecule(const Molecule& m);
  Molecule(std::istream& istrm, streamtype type, const std::vector<Component>* components=NULL);

  ~Molecule()
    { assert(m_sites_d); delete[] m_sites_d; assert(m_osites_e); delete[] m_osites_e; assert(m_sites_F); delete[] m_sites_F; }
  /** get the ID */
  unsigned long id() const { return m_id; }
  void setid(unsigned long id) { this->m_id = id; }
  /** get the Component */
  int componentid() const { return m_componentid; }
  /** get the position */
  double r(unsigned short d) const { return m_r[d]; }

  double oldr (unsigned short d) const { return m_oldr[d];}

  /** get the velocity */
  double v(unsigned short d) const { return m_v[d]; }
  /** get the Orientation */
  const Quaternion& q() const { return m_q; }

  inline void move(int d, double dr) { m_r[d] += dr; }
  
  /** get the rotatational speed */
  double D(unsigned short d) const { return m_D[d]; }

  /** get F */
  double F(unsigned short d) const {return m_F[d]; }
  /** get M */
  double M(unsigned short d) const {return m_M[d]; }

  //double Upot() const { return m_Upot; }
  double Utrans() const { return .5*m_m*(m_v[0]*m_v[0]+m_v[1]*m_v[1]+m_v[2]*m_v[2]); }
  double Urot();

  /** get number of sites */
  unsigned int numSites() const { return m_numsites; }
  unsigned int numLJcenters() const { return m_ljcenters->size(); }
  unsigned int numCharges() const { return m_charges->size(); }
  unsigned int numDipoles() const { return m_dipoles->size(); }
  unsigned int numQuadrupoles() const { return m_quadrupoles->size(); }
  unsigned int numTersoff() const { assert(m_tersoff); return m_tersoff->size(); }

  const double* site_d(unsigned int i) const { return &(m_sites_d[3*i]); }
  const double* osite_e(unsigned int i) const { return &(m_osites_e[3*i]); }
  const double* site_F(unsigned int i) const { return &(m_sites_F[3*i]); }
  const double* ljcenter_d(unsigned int i) const { return &(m_ljcenters_d[3*i]); }
  const double* ljcenter_F(unsigned int i) const { return &(m_ljcenters_F[3*i]); }
  const double* charge_d(unsigned int i) const { return &(m_charges_d[3*i]); }
  const double* charge_F(unsigned int i) const { return &(m_charges_F[3*i]); }
  const double* dipole_d(unsigned int i) const { return &(m_dipoles_d[3*i]); }
  const double* dipole_e(unsigned int i) const { return &(m_dipoles_e[3*i]); }
  const double* dipole_F(unsigned int i) const { return &(m_dipoles_F[3*i]); }
  const double* quadrupole_d(unsigned int i) const { return &(m_quadrupoles_d[3*i]); }
  const double* quadrupole_e(unsigned int i) const { return &(m_quadrupoles_e[3*i]); }
  const double* quadrupole_F(unsigned int i) const { return &(m_quadrupoles_F[3*i]); }
  const double* tersoff_d(unsigned int i) const { return &(m_tersoff_d[3*i]); }
  const double* tersoff_F(unsigned int i) const { return &(m_tersoff_F[3*i]); }

  /** get object memory size */
  static unsigned long memsize() { return sizeof(Molecule); }

  /** set the position */
  void setr(unsigned short d, double r) { m_r[d]=r; }

  void setOldRFromR () {for (int i = 0; i < 3; i++) m_oldr[i] = m_r[i];}

  /** calculate the difference vector and return the square (euclidean) distance */
  double dist2(const Molecule& a, double dr[]) const
    { double d2=0.; for(unsigned short d=0;d<3;++d) { dr[d]=a.m_r[d]-m_r[d]; d2+=dr[d]*dr[d]; } return d2; }
  double dist2(const Molecule& a, double L[3], double dr[]) const;
  /** calculate and return the square velocity */
  double v2() const {return m_v[0]*m_v[0]+m_v[1]*m_v[1]+m_v[2]*m_v[2]; }

  void setFM(double Fx, double Fy, double Fz, double Mx, double My, double Mz)
    { m_F[0]=Fx; m_F[1]=Fy; m_F[2]=Fz; m_M[0]=Mx; m_M[1]=My; m_M[2]=My; }
  void scale_v(double s) { for(unsigned short d=0;d<3;++d) m_v[d]*=s; }
  void scale_v(double s, double offx, double offy, double offz);
  void scale_F(double s) { for(unsigned short d=0;d<3;++d) m_F[d]*=s; }
  void scale_D(double s) { for(unsigned short d=0;d<3;++d) m_D[d]*=s; }
  void scale_M(double s) { for(unsigned short d=0;d<3;++d) m_M[d]*=s; }

  void Fadd(const double a[]) { for(unsigned short d=0;d<3;++d) m_F[d]+=a[d]; }
  void Fsub(const double a[]) { for(unsigned short d=0;d<3;++d) m_F[d]-=a[d]; }

  void Madd(const double a[]) { for(unsigned short d=0;d<3;++d) m_M[d]+=a[d]; }
  void Msub(const double a[]) { for(unsigned short d=0;d<3;++d) m_M[d]-=a[d]; }

  void vadd(const double ax, const double ay, const double az)
  {
     m_v[0] += ax; m_v[1] += ay; m_v[2] += az;
  }
  void vsub(const double ax, const double ay, const double az)
  {
     m_v[0] -= ax; m_v[1] -= ay; m_v[2] -= az;
  }
  void setXY() { fixedx = m_r[0]; fixedy = m_r[1]; }
  void resetXY()
  {
     m_v[0] = 0.0;
     m_v[1] = 0.0;
     m_F[1] = 0.0;
     m_F[0] = 0.0; 
     m_r[0] = fixedx;
     m_r[1] = fixedy;
  }
  
  void Fsiteadd(unsigned int i, double a[])
    { double* Fsite=&(m_sites_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fsitesub(unsigned int i, double a[])
    { double* Fsite=&(m_sites_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fljcenteradd(unsigned int i, double a[])
    { double* Fsite=&(m_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fljcentersub(unsigned int i, double a[])
    { double* Fsite=&(m_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fchargeadd(unsigned int i, double a[])
    { double* Fsite=&(m_charges_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fchargesub(unsigned int i, double a[])
    { double* Fsite=&(m_charges_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fdipoleadd(unsigned int i, double a[])
    { double* Fsite=&(m_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fdipolesub(unsigned int i, double a[])
    { double* Fsite=&(m_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fquadrupoleadd(unsigned int i, double a[])
    { double* Fsite=&(m_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fquadrupolesub(unsigned int i, double a[])
    { double* Fsite=&(m_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Ftersoffadd(unsigned int i, double a[])
    { double* Fsite=&(m_tersoff_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Ftersoffsub(unsigned int i, double a[])
    { double* Fsite=&(m_tersoff_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }

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
  
  inline unsigned getCurTN() { return this->m_curTN; }
  inline Molecule* getTersoffNeighbour(unsigned i) { return this->m_Tersoff_neighbours_first[i]; }
  inline bool getPairCode(unsigned i) { return this->m_Tersoff_neighbours_second[i]; }
  inline void clearTersoffNeighbourList() { this->m_curTN = 0; }
  void addTersoffNeighbour(Molecule* m, bool pairType);
  double tersoffParameters(double params[15]); //returns delta_r

  // clear forces and moments
  void clearFM();
  void check(unsigned long id);
  
private:
  
  static Domain* _domain;
  unsigned long m_id; // IDentification number of that molecule
  int m_componentid;  // IDentification number of its component type
  double m_r[3];  // position coordinates
  double m_oldr[3]; // position coordinates last step
  double m_v[3];  // velocity
  Quaternion m_q; // orientation
  double m_D[3];  // angular momentum

  double m_F[3];  // forces
  double m_M[3];  // moments

  const std::vector<LJcenter>* m_ljcenters;
  const std::vector<Charge>* m_charges;
  const std::vector<Dipole>* m_dipoles;
  const std::vector<Quadrupole>* m_quadrupoles;
  const std::vector<Tersoff>* m_tersoff;
  
  double m_m; // total mass
  double m_I[3],m_invI[3];  // moment of inertia for principal axes and it's inverse
  std::size_t m_numsites; // number of sites
  std::size_t m_numorientedsites; // number of oriented sites (subset of sites)
  // global site coordinates relative to site origin
  // row order: dx1,dy1,dz1,dx2,dy2,dz2,...
  double *m_sites_d;
  double *m_ljcenters_d, *m_charges_d, *m_dipoles_d,
         *m_quadrupoles_d, *m_tersoff_d;
  // site orientation
  double *m_osites_e;
  double *m_dipoles_e, *m_quadrupoles_e;
  // site Forces
  // row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
  double* m_sites_F;
  double *m_ljcenters_F, *m_charges_F, *m_dipoles_F,
         *m_quadrupoles_F, *m_tersoff_F;

  Molecule* m_Tersoff_neighbours_first[MAXTN];
  bool m_Tersoff_neighbours_second[MAXTN];
  int m_curTN;
  double fixedx, fixedy;

  // setup cache values/properties
  void setupCache(const std::vector<Component>* components);
  // calculate forces and moments for already given site forces
  void calcFM();
};

#endif /*MOLECULE_H_*/
