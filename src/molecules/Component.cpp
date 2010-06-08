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

#include "molecules/Component.h"
#include <iostream>
#include <iomanip>

using namespace std;

Component::Component(unsigned int id)
{
  m_id=id;
  m_m=0.;
  m_I[0]=m_I[1]=m_I[2]=m_I[3]=m_I[4]=m_I[5]=0.;
  m_rot_dof=0;
  m_Ipa[0]=m_Ipa[1]=m_Ipa[2]=0.;
  m_numMolecules=0;
  this->maximalTersoffExternalRadius = 0.0;

  m_ljcenters = vector<LJcenter>();
  m_charges = vector<Charge>();
  m_quadrupoles = vector<Quadrupole>();
  m_dipoles = vector<Dipole>();
  m_tersoff = vector<Tersoff>();
}

void Component::setI(double Ixx,double Iyy,double Izz
                    ,double Ixy,double Ixz,double Iyz)
{
  m_I[0]=Ixx;
  m_I[1]=Iyy;
  m_I[2]=Izz;
  m_I[3]=Ixy;
  m_I[4]=Ixz;
  m_I[5]=Iyz;
}

void Component::addI(double Ixx,double Iyy,double Izz
                    ,double Ixy,double Ixz,double Iyz)
{
  m_I[0]+=Ixx;
  m_I[1]+=Iyy;
  m_I[2]+=Izz;
  m_I[3]+=Ixy;
  m_I[4]+=Ixz;
  m_I[5]+=Iyz;
}

void Component::addLJcenter(
   double x, double y, double z, double m, double eps, double sigma,
   double rc, bool TRUNCATED_SHIFTED
) {
  if(TRUNCATED_SHIFTED)
  {
     double sigperrc2 = sigma*sigma/(rc*rc);
     double sigperrc6 = sigperrc2*sigperrc2*sigperrc2;
     double shift6 = 24.0*eps * (sigperrc6 - sigperrc6*sigperrc6);
     m_ljcenters.push_back(LJcenter(x, y, z, m, eps, sigma, rc, shift6));
  }
  else
  {
     m_ljcenters.push_back(LJcenter(x, y, z, m, eps, sigma, rc, 0.0)); 
  }
  
  m_m+=m;
  // assume the input is already transformed to the principal axes system
  // (and therefore the origin is the center of mass)
  m_I[0]+=m*(y*y+z*z);
  m_I[1]+=m*(x*x+z*z);
  m_I[2]+=m*(x*x+y*y);
  m_I[3]-=m*x*y;
  m_I[4]-=m*x*z;
  m_I[5]-=m*y*z;
  
  m_rot_dof=3;
  for(unsigned short d=0;d<3;++d)
  {
    m_Ipa[d]=m_I[d];
    if(m_Ipa[d]==0.) --m_rot_dof;
  }
}

void Component::addCharge(double x, double y, double z, double m, double q)
{
   m_charges.push_back( Charge(x, y, z, m, q) );
   m_m += m;

   // assume the input is already transformed to the principal axes system
   // (and therefore the origin is the center of mass)
   m_I[0]+=m*(y*y+z*z);
   m_I[1]+=m*(x*x+z*z);
   m_I[2]+=m*(x*x+y*y);
   m_I[3]-=m*x*y;
   m_I[4]-=m*x*z;
   m_I[5]-=m*y*z;

   m_rot_dof=3;
   for(unsigned short d=0;d<3;++d)
   {
      m_Ipa[d]=m_I[d];
      if(m_Ipa[d]==0.) --m_rot_dof;
   }
}

void Component::addDipole(double x, double y, double z
                         ,double eMyx, double eMyy, double eMyz, double eMyabs)
{
  m_dipoles.push_back(Dipole(x,y,z,eMyx,eMyy,eMyz,eMyabs));
  // massless...
}

void Component::addQuadrupole(double x, double y, double z
                             ,double eQx, double eQy, double eQz, double eQabs)
{
  m_quadrupoles.push_back(Quadrupole(x,y,z,eQx,eQy,eQz,eQabs));
  // massless...
}

void Component::addTersoff( double x, double y, double z,
                            double m, double A, double B, double lambda, double mu, double R,
                            double S, double c, double d, double h, double n, double beta  )
{
   if(S > this->maximalTersoffExternalRadius) maximalTersoffExternalRadius = S;
   m_tersoff.push_back( Tersoff(x, y, z, m, A, B, lambda, mu, R, S, c, d, h, n, beta) );

   m_m+=m;
   // assume the input is already transformed to the principal axes system
   // (and therefore the origin is the center of mass)
   m_I[0]+=m*(y*y+z*z);
   m_I[1]+=m*(x*x+z*z);
   m_I[2]+=m*(x*x+y*y);
   m_I[3]-=m*x*y;
   m_I[4]-=m*x*z;
   m_I[5]-=m*y*z;

   m_rot_dof=3;
   for(unsigned short d=0;d<3;++d)
   {
     m_Ipa[d]=m_I[d];
     if(m_Ipa[d]==0.) --m_rot_dof;
   }
}

void Component::write(std::ostream& ostrm) const
{
  ostrm << m_ljcenters.size() << "\t" << m_charges.size() << "\t"
        << m_dipoles.size() << "\t" << m_quadrupoles.size() << "\t"
	<< m_tersoff.size() << "\n";
  for(std::vector<LJcenter>::const_iterator pos=m_ljcenters.begin();pos!=m_ljcenters.end();++pos)
  {
    pos->write(ostrm);
    ostrm << endl;
  }
  for(std::vector<Charge>::const_iterator pos=m_charges.begin();pos!=m_charges.end();++pos)
  {
    pos->write(ostrm);
    ostrm << endl;
  }
  for(std::vector<Dipole>::const_iterator pos=m_dipoles.begin();pos!=m_dipoles.end();++pos)
  {
    pos->write(ostrm);
    ostrm << endl;
  }
  for(std::vector<Quadrupole>::const_iterator pos=m_quadrupoles.begin();pos!=m_quadrupoles.end();++pos)
  {
    pos->write(ostrm);
    ostrm << endl;
  }
  for(std::vector<Tersoff>::const_iterator pos=m_tersoff.begin();pos!=m_tersoff.end();++pos)
  {
    pos->write(ostrm);
    ostrm << endl;
  }
  ostrm << m_Ipa[0] << " " << m_Ipa[1] << " " << m_Ipa[2] << endl;
}

void Component::writePOVobjs(std::ostream& ostrm, string para) const
{
  if(numLJcenters()<=0) return;
  if(numLJcenters()==1)
  {
    ostrm << "sphere {<" << m_ljcenters.front().rx() << "," << m_ljcenters.front().ry() << "," << m_ljcenters.front().rz() << ">," << .5*m_ljcenters.front().sigma() << " " << para << "}";
  }
  else
  {
    ostrm << "blob { threshold 0.01 ";
    for(std::vector<LJcenter>::const_iterator pos=m_ljcenters.begin();pos!=m_ljcenters.end();++pos)
      ostrm << "sphere {<" << pos->rx() << "," << pos->ry() << "," << pos->rz() << ">," << .5*pos->sigma() << ", strength 1 } ";
    ostrm << para << "}";
  }
  ostrm << flush;
}

void Component::writeVIM(std::ostream& ostrm)
{
   for(std::vector<LJcenter>::const_iterator pos=m_ljcenters.begin();pos!=m_ljcenters.end();++pos)
   {
      ostrm << "~ " << this->m_id+1 << " LJ "
            << setw(7) << pos->rx() << ' ' << setw(7) << pos->ry() << ' ' << setw(7) << pos->rz() << ' '
            << setw(6) << pos->sigma() << ' ' << setw(2) << (1 + (this->m_id % 9)) << "\n";
   }
   ostrm << flush;
}
