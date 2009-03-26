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
 ****************c***********************************************************/
#ifndef MOLECULE_H_
#define MOLECULE_H_

/**
Molecule modeled as LJ sphere

@author Martin Bernreuther
*/

#include "molecules/Quaternion.h"
#include "molecules/Comp2Param.h"
#include "molecules/Site.h"
#include "molecules/Component.h"
#include "integrators/Integrator.h"
class Domain;

#include <vector>
#include <iostream>
#include <string>
using namespace std;

#include <cassert>

class Molecule{
//friend class integrators::Integrator;
//friend class integrators::Leapfrog;
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
  /** get the Component */
  int componentid() const { return m_componentid; }
  /** get the position */
  double r(unsigned short d) const { return m_r[d]; }
  /** get the velocity */
  double v(unsigned short d) const { return m_v[d]; }
  /** get the Orientation */
  const Quaternion& q() const { return m_q; }

  /** get the rotatational speed */
  double D(unsigned short d) const { return m_D[d]; }

  /** get F */
  double F(unsigned short d) const {return m_F[d]; }
  /** get M */
  double M(unsigned short d) const {return m_M[d]; }

  //double Upot() const { return m_Upot; }
  double Ukin() const { return .5*m_m*(m_v[0]*m_v[0]+m_v[1]*m_v[1]+m_v[2]*m_v[2]); }

  /** get number of sites */
  unsigned int numSites() const { return m_numsites; }
  unsigned int numLJcenters() const { return m_ljcenters->size(); }
  unsigned int numDipoles() const { return m_dipoles->size(); }
  unsigned int numQuadrupoles() const { return m_quadrupoles->size(); }

  const double* site_d(unsigned int i) const { return &(m_sites_d[3*i]); }
  const double* osite_e(unsigned int i) const { return &(m_osites_e[3*i]); }
  const double* site_F(unsigned int i) const { return &(m_sites_F[3*i]); }
  const double* ljcenter_d(unsigned int i) const { return &(m_ljcenters_d[3*i]); }
  const double* ljcenter_F(unsigned int i) const { return &(m_ljcenters_F[3*i]); }
  const double* dipole_d(unsigned int i) const { return &(m_dipoles_d[3*i]); }
  const double* dipole_e(unsigned int i) const { return &(m_dipoles_e[3*i]); }
  const double* dipole_F(unsigned int i) const { return &(m_dipoles_F[3*i]); }
  const double* quadrupole_d(unsigned int i) const { return &(m_quadrupoles_d[3*i]); }
  const double* quadrupole_e(unsigned int i) const { return &(m_quadrupoles_e[3*i]); }
  const double* quadrupole_F(unsigned int i) const { return &(m_quadrupoles_F[3*i]); }

  /** get object memory size */
  static unsigned long memsize() { return sizeof(Molecule); }

  /** set the position */
  void setr(unsigned short d, double r) { m_r[d]=r; }

  /** calculate the difference vector and return the square (euclidean) distance */
  double dist2(const Molecule& a, double dr[]) const
    { double d2=0.; for(unsigned short d=0;d<3;++d) { dr[d]=a.m_r[d]-m_r[d]; d2+=dr[d]*dr[d]; } return d2; }
  double dist2(const Molecule& a, double L[3], double dr[]) const;
  /** calculate and return the square velocity */
  double v2() const {return m_v[0]*m_v[0]+m_v[1]*m_v[1]+m_v[2]*m_v[2]; }

  void setFM(double Fx, double Fy, double Fz, double Mx, double My, double Mz)
    { m_F[0]=Fx; m_F[1]=Fy; m_F[2]=Fz; m_M[0]=Mx; m_M[1]=My; m_M[2]=My; }
  void scale_v(double s) { for(unsigned short d=0;d<3;++d) m_v[d]*=s; }
  void scale_F(double s) { for(unsigned short d=0;d<3;++d) m_F[d]*=s; }
  void scale_D(double s) { for(unsigned short d=0;d<3;++d) m_D[d]*=s; }

  void Fadd(const double a[]) { for(unsigned short d=0;d<3;++d) m_F[d]+=a[d]; }
  void Fsub(const double a[]) { for(unsigned short d=0;d<3;++d) m_F[d]-=a[d]; }

  void Madd(const double a[]) { for(unsigned short d=0;d<3;++d) m_M[d]+=a[d]; }
  void Msub(const double a[]) { for(unsigned short d=0;d<3;++d) m_M[d]-=a[d]; }

  void Fsiteadd(unsigned int i, double a[])
    { double* Fsite=&(m_sites_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fsitesub(unsigned int i, double a[])
    { double* Fsite=&(m_sites_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fljcenteradd(unsigned int i, double a[])
    { double* Fsite=&(m_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fljcentersub(unsigned int i, double a[])
    { double* Fsite=&(m_ljcenters_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fdipoleadd(unsigned int i, double a[])
    { double* Fsite=&(m_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fdipolesub(unsigned int i, double a[])
    { double* Fsite=&(m_dipoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }
  void Fquadrupoleadd(unsigned int i, double a[])
    { double* Fsite=&(m_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]+=a[d]; }
  void Fquadrupolesub(unsigned int i, double a[])
    { double* Fsite=&(m_quadrupoles_F[3*i]); for(unsigned short d=0;d<3;++d) Fsite[d]-=a[d]; }

  //void Upotadd(double u) { m_Upot+=u; }
  //void Upotsub(double u) { m_Upot-=u; }

  void upd_preF(double dt, double vcorr=1., double Dcorr=1.);
  void upd_cache();
  void upd_postF(double dt_halve, double& summv2, double& sumIw2);

  //! calculate summv2 and sumIw2
  //! @todo what is sumIw2?
  //! @todo comment
  void calculate_mv2_Iw2(double& summv2, double& sumIw2);

  /** write information to stream */
  void write(std::ostream& ostrm) const;
  /** write binary information to stream */
  void save_restart(std::ostream& ostrm) const;

//  Linked cells
  Molecule* nextinCell() const { return m_nextincell; }
  
  static void setDomain(Domain* domain);
  static void setblubb(double a);


private:
  
  static Domain* _domain;
  unsigned long m_id; // IDentification number of that molecule
  int m_componentid;  // IDentification number of its component type
  double m_r[3];  // position coordinates
  double m_v[3];  // velocity
  Quaternion m_q; // orientation
  double m_D[3];  // angular momentum

  double m_F[3];  // forces
  double m_M[3];  // moments
  //double m_Upot;  // potential energy

//  caching component values
  //const std::vector<Site>* m_sites; // doesn't work that way... see component.h
  const std::vector<LJcenter>* m_ljcenters;
  const std::vector<Dipole>* m_dipoles;
  const std::vector<Quadrupole>* m_quadrupoles;
  double m_m; // total mass
  double m_I[3],m_invI[3];  // moment of inertia for principal axes and it's inverse
  std::size_t m_numsites; // number of sites
  std::size_t m_numorientedsites; // number of oriented sites (subset of sites)
  // global site coordinates relative to site origin
  // row order: dx1,dy1,dz1,dx2,dy2,dz2,...
  double *m_sites_d;
  double *m_ljcenters_d, *m_dipoles_d, *m_quadrupoles_d;
  // site orientation
  double *m_osites_e;
  double *m_dipoles_e, *m_quadrupoles_e;
  // site Forces
  // row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
  double* m_sites_F;
  double *m_ljcenters_F, *m_dipoles_F, *m_quadrupoles_F;

  // used by CELL data structure for intrusive single-linked "collision" lists
  Molecule* m_nextincell;

  // used by Cells::addMolecule to modify m_nextincell
  void setNextinCell(Molecule* next) { m_nextincell=next; }
//  friend void Cells::addMolecule(Molecule* atom);


  
  // setup cache values/properties
  void setupCache(const std::vector<Component>* components);
  // clear forces and moments
  void clearFM();
  // calculate forces and moments for already given site forces
  void calcFM();
  

};

#endif /*MOLECULE_H_*/
