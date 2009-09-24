/***************************************************************************
 *   Copyright (C) 2009 by Martin Bernreuther et al.                       *
 *   bernreuther@hlrs.de                                                   *
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
#include <math.h>

/** Site
@author Martin Bernreuther et al. (2009)
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
@author Martin Bernreuther et al.
*/
class LJcenter : public Site
{
public:
  /// Constructor: pass shift = 0.0 for the full LJ potential
  LJcenter(double x, double y, double z ,double m, double eps, double sigma, double rc, double shift)
    : Site(x,y,z,m), m_eps(eps), m_sigma(sigma), m_rc(rc)
    { m_r[0]=x; m_r[1]=y; m_r[2]=z; uLJshift6 = shift; }
  /// Constructor reading from stream
  LJcenter(std::istream& istrm) { istrm >> m_r[0] >> m_r[1] >> m_r[2] >> m_m >> m_eps >> m_sigma >> uLJshift6; }
  /// write to stream
  void write(std::ostream& ostrm) const { ostrm << m_r[0] << " " << m_r[1] << " "  << m_r[2]<< "\t"  << m_m << "\t"  << m_eps << " " << m_sigma << " " << m_rc << " " << uLJshift6; }
  /// get strength
  double eps() const { return m_eps; }
  /// get diameter
  double sigma() const { return m_sigma; }
  /// get truncation option
  bool TRUNCATED_SHIFTED() { return (this->uLJshift6 != 0.0); }
  double shift6() const { return this->uLJshift6; }

private:
  double m_eps; // strength
  double m_sigma; // diameter
  double m_rc; // cutoff radius
  double uLJshift6;  // truncation offset of the LJ potential
};

class Charge: public Site
{
 public:
   /// Constructor
   Charge(double x, double y, double z, double m, double q)
      : Site(x, y, z, m), m_q(q)
   {
      m_r[0]=x; m_r[1]=y; m_r[2]=z;
      m_m = m;
      m_q = q;
   }

   /// Constructor reading from stream
   Charge(std::istream& istrm)
   {
      istrm >> m_r[0] >> m_r[1] >> m_r[2] >> m_m >> m_q;
   }
   /// write to stream
   void write(std::ostream& ostrm) const 
   {
      ostrm << m_r[0] << " " << m_r[1] << " " << m_r[2] << "\t" << m_m << " " << m_q;
   }
   /// get charge
   double q() const { return m_q; }

 private:
   /// charge
   double m_q;
};

class Tersoff: public Site
{
 public:
   Tersoff( double x, double y, double z,
            double m, double A, double B,
            double lambda, double mu, double R, double S,
            double c, double d, double h, double n, double beta )
      : Site(x, y, z, m)
   {
      this->m_r[0] = x;
      this->m_r[1] = y;
      this->m_r[2] = z;

      this->m_m = m;
      this->m_A = A;
      this->m_B = B;

      this->m_minus_lambda = -lambda;
      this->m_minus_mu = -mu;
      this->m_R = R;
      this->m_S = S;

      this->m_c_square = c*c;
      this->m_d_square = d*d;
      this->m_h = h;
      this->m_n = n;
      this->m_beta = beta;
   }

   /// Constructor reading from stream
   Tersoff(std::istream& istrm)
   {
      double lambda, mu, c, d;
      istrm >> m_r[0] >> m_r[1] >> m_r[2]
            >> m_m >> m_A >> m_B >> lambda >> mu
            >> m_R >> m_S >> c
            >> d >> m_h >> m_n >> m_beta;
      m_minus_lambda = -lambda;
      m_minus_mu = -mu;
      m_c_square = c*c;
      m_d_square = d*d;
   }

   //! @brief write to stream
   //!
   void write(std::ostream& ostrm) const
   {
      ostrm << m_r[0] << " " << m_r[1] << " " << m_r[2] << "\t"
            << m_m << "\t"  << m_A << " "  << m_B << " " << -m_minus_lambda << " " << -m_minus_mu << "\t"
            << m_R << " " << m_S << "\t" << sqrt(m_c_square) << " "
            << sqrt(m_d_square) << "\t" << m_h << " " << m_n << " " << m_beta;
   }

   //! @brief get repulsive coefficient
   //!
   double A() const { return this->m_A; }

   //! @brief get attractive coefficient
   //!
   double B() const { return this->m_B; }

   //! @brief get repulsive coexponent
   //!
   double minusLambda() const { return this->m_minus_lambda; }

   //! @brief get attractive coexponent
   //!
   double minusMu() const { return this->m_minus_mu; }

   //! @brief get internal radius
   //!
   double R() const { return this->m_R; }

   //! @brief get external radius
   //!
   double S() const { return this->m_S; }

   //! @brief get c square Tersoff interaction parameter
   //!
   double cSquare() const { return this->m_c_square; }

   //! @brief get d square Tersoff interaction parameter
   //!
   double dSquare() const { return this->m_d_square; }

   //! @brief get h Tersoff interaction parameter
   //!
   double h() const { return this->m_h; }

   //! @brief get n Tersoff interaction parameter
   //!
   double n() const { return this->m_n; }

   //! @brief get beta Tersoff interaction parameter
   //!
   double beta() const { return this->m_beta; }

 private:
   //! repulsive coefficient
   double m_A;
   //! attractive coefficient
   double m_B;
   //! repulsive coexponent
   double m_minus_lambda;
   //! attractive coexponent
   double m_minus_mu;

   //! internal radius
   double m_R;
   //! external radius
   double m_S;

   //! c square Tersoff interaction parameter
   double m_c_square;
   //! d square Tersoff interaction parameter
   double m_d_square;
   //! h Tersoff interaction parameter
   double m_h;
   //! n Tersoff interaction parameter
   double m_n;
   //! beta Tersoff interaction parameter
   double m_beta;
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
