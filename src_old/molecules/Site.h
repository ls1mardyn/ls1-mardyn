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
#ifndef SITE_H_
#define SITE_H_

#include <iostream>

/** Site
@author Martin Bernreuther
*/
class Site
{
public:
  /// get x-coordinate of position
  double rx() const { return m_r[0]; }
  /// get y-coordinate of position
  double ry() const { return m_r[1]; }
  /// get z-coordinate of position
  double rz() const { return m_r[2]; }
  /// get position vector
  const double* r() const { return m_r; }
  /// get mass
  double m() const { return m_m; }
  /// translate coordinates to new origin
  void translateOrigin(double neworigin[3])
    { for(unsigned short d=0;d<3;++d) m_r[d]-=neworigin[d]; }

protected:
  /// Constructor
  Site(double x=0., double y=0., double z=0., double m=0.)
    : m_m(m)
    { m_r[0]=x; m_r[1]=y; m_r[2]=z; }

//private:
  double m_r[3];  // position coordinates
  double m_m; // mass
};


/** LJcenter
    Lennard-Jones 12-6 center
@author Martin Bernreuther
*/
class LJcenter : public Site
{
public:
  /// Constructor
  LJcenter(double x, double y, double z ,double m, double eps, double sigma)
    : Site(x,y,z,m), m_eps(eps), m_sigma(sigma)
    { m_r[0]=x; m_r[1]=y; m_r[2]=z; }
  /// Constructor reading from stream
  LJcenter(std::istream& istrm) { istrm >> m_r[0] >> m_r[1] >> m_r[2] >> m_m >> m_eps >> m_sigma; }
  /// write to stream
  void write(std::ostream& ostrm) const { ostrm << m_r[0] << " " << m_r[1] << " "  << m_r[2]<< "\t"  << m_m << "\t"  << m_eps << " "  << m_sigma; }
  /// get strength
  double eps() const { return m_eps; }
  /// get diameter
  double sigma() const { return m_sigma; }

private:
  double m_eps; // strength
  double m_sigma; // diameter
};


/** OrientedSite
@author Martin Bernreuther
*/
class OrientedSite : public Site
{
public:
  double ex() const { return m_e[0]; }
  double ey() const { return m_e[1]; }
  double ez() const { return m_e[2]; }
  const double* e() const { return m_e; }

protected:
  /// Constructor
  OrientedSite(double x=0., double y=0., double z=0., double m=0., double ex=0., double ey=0., double ez=0.)
        : Site(x,y,z,m)
    { m_e[0]=ex; m_e[1]=ey; m_e[2]=ez; }

//private:
  double m_e[3];
};


/** Dipole
@author Martin Bernreuther
*/
class Dipole : public OrientedSite
{
public:
  /// Constructor
  Dipole(double x, double y, double z , double eMyx, double eMyy, double eMyz, double absMy)
    : OrientedSite(x,y,z,0.,eMyx,eMyy,eMyz), m_absMy(absMy) { }
  /// Constructor reading from stream
  Dipole(std::istream& istrm) { istrm >> m_r[0] >> m_r[1] >> m_r[2] >> m_e[0] >> m_e[1] >> m_e[2] >> m_absMy; m_m=0.;}
  /// write to stream
  void write(std::ostream& ostrm) const { ostrm << m_r[0] << " " << m_r[1] << " " << m_r[2] << "\t" << m_e[0] << " " << m_e[1] << " " << m_e[2] << "\t" << m_absMy; }
  double absMy() const { return m_absMy; }

private:
  double m_absMy;
};

/** Quadrupole
@author Martin Bernreuther
*/
class Quadrupole : public OrientedSite
{
public:
  /// Constructor
  Quadrupole(double x, double y, double z , double eQx, double eQy, double eQz, double absQ)
      : OrientedSite(x,y,z,0.,eQx,eQy,eQz), m_absQ(absQ) { }
  /// Constructor reading from stream
  Quadrupole(std::istream& istrm) { istrm >> m_r[0] >> m_r[1] >> m_r[2] >> m_e[0] >> m_e[1] >> m_e[2] >> m_absQ; m_m=0.;}
  /// write to stream
  void write(std::ostream& ostrm) const { ostrm << m_r[0] << " " << m_r[1] << " " << m_r[2] << "\t" << m_e[0] << " " << m_e[1] << " " << m_e[2] << " " << m_absQ; }
  double absQ() const { return m_absQ; }

private:
  double m_absQ;
};

#endif /*SITE_H_*/
