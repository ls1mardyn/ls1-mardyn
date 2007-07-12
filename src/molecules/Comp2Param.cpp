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

#include <cmath>

#include "molecules/Comp2Param.h"

using namespace std;

void Comp2Param::initialize(const vector<Component>& components, const vector<double>& mixcoeff
                           , double epsRF, double rc)
{
  m_numcomp=components.size();
  m_ssparatbl.redim(m_numcomp,m_numcomp);

  // interaction between LJ centers
  vector<double>::const_iterator mixpos=mixcoeff.begin();
  for(unsigned int compi=0;compi<m_numcomp;++compi)
  {
    ParaStrm& pstrmii=m_ssparatbl(compi,compi);
    unsigned int nci=components[compi].numLJcenters();
    double epsi,sigi,epsj,sigj,epsilon24,sigma2;
    // interaction between same components
    for(unsigned int centeri=0;centeri<nci;++centeri)
    {
      const LJcenter& ljcenteri=static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
      epsi=ljcenteri.eps();
      sigi=ljcenteri.sigma();
      for(unsigned int centerj=0;centerj<nci;++centerj)
      {
        const LJcenter& ljcenterj=static_cast<const LJcenter&>(components[compi].ljcenter(centerj));
        epsj=ljcenterj.eps();
        sigj=ljcenterj.sigma();
        epsilon24=24.*sqrt(epsi*epsj);
        sigma2=.5*(sigi+sigj); sigma2*=sigma2;
        pstrmii << epsilon24;
        pstrmii << sigma2;
      }
    }
    // interaction between different components
    for(unsigned int compj=compi+1;compj<m_numcomp;++compj)
    {
      ParaStrm& pstrmij=m_ssparatbl(compi,compj);
      unsigned int ncj=components[compj].numLJcenters();
      double xi=*mixpos; ++mixpos;
      double eta=*mixpos; ++mixpos;
      for(unsigned int centeri=0;centeri<nci;++centeri)
      {
        const LJcenter& ljcenteri=static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
        epsi=ljcenteri.eps();
        sigi=ljcenteri.sigma();
        for(unsigned int centerj=0;centerj<ncj;++centerj)
        {
          const LJcenter& ljcenterj=static_cast<const LJcenter&>(components[compj].ljcenter(centerj));
          epsj=ljcenterj.eps();
          sigj=ljcenterj.sigma();
          epsilon24=24.*xi*sqrt(epsi*epsj);
          sigma2=eta*.5*(sigi+sigj); sigma2*=sigma2;
          pstrmij << epsilon24;
          pstrmij << sigma2;
        }
      }
      ParaStrm& pstrmji=m_ssparatbl(compj,compi);
      for(unsigned int centerj=0;centerj<ncj;++centerj)
      {
        const LJcenter& ljcenterj=static_cast<const LJcenter&>(components[compj].ljcenter(centerj));
        epsj=ljcenterj.eps();
        sigj=ljcenterj.sigma();
        for(unsigned int centeri=0;centeri<nci;++centeri)
        {
          const LJcenter& ljcenteri=static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
          epsi=ljcenteri.eps();
          sigi=ljcenteri.sigma();
          epsilon24=24.*xi*sqrt(epsi*epsj);
          sigma2=eta*.5*(sigi+sigj); sigma2*=sigma2;
          pstrmji << epsilon24;
          pstrmji << sigma2;
        }
      }
    }
  }

  for(unsigned int compi=0;compi<m_numcomp;++compi)
  {
    unsigned int ndi=components[compi].numDipoles();
    unsigned int nqi=components[compi].numQuadrupoles();
    for(unsigned int compj=0;compj<m_numcomp;++compj)
    {
      ParaStrm& pstrmij=m_ssparatbl(compi,compj);
      unsigned int ndj=components[compj].numDipoles();
      unsigned int nqj=components[compj].numQuadrupoles();
      for(unsigned int di=0;di<ndi;++di)
      {
        const Dipole& dipolei=static_cast<const Dipole&>(components[compi].dipole(di));
        double absmyi=dipolei.absMy();
        double epsRFInvrc3=2.*(epsRF-1.)/((rc*rc*rc)*(2.*epsRF+1.));
        // Dipole-Dipole
        for(unsigned int dj=0;dj<ndj;++dj)
        {
          const Dipole& dipolej=static_cast<const Dipole&>(components[compj].dipole(dj));
          double absmyj=dipolej.absMy();
          double my2=absmyi*absmyj;
          pstrmij << my2;
          double rffac=my2*epsRFInvrc3;
          pstrmij << rffac;
        }
        // Dipole-Quadrupole
        for(unsigned int qj=0;qj<nqj;++qj)
        {
          const Quadrupole& quadrupolej=static_cast<const Quadrupole&>(components[compj].quadrupole(qj));
          double absqj=quadrupolej.absQ();
          double myq15=1.5*absmyi*absqj;
          pstrmij << myq15;
        }
      }
      for(unsigned int qi=0;qi<nqi;++qi)
      {
        const Quadrupole& quadrupolei=static_cast<const Quadrupole&>(components[compi].quadrupole(qi));
        double absqi=quadrupolei.absQ();
        // Quadrupole-Dipole
        for(unsigned int dj=0;dj<ndj;++dj)
        {
          const Dipole& dipolej=static_cast<const Dipole&>(components[compj].dipole(dj));
          double absmyj=dipolej.absMy();
          double qmy15=1.5*absqi*absmyj;
          pstrmij << qmy15;
        }
        // Quadrupole-Quadrupole
        for(unsigned int qj=0;qj<nqj;++qj)
        {
          const Quadrupole& quadrupolej=static_cast<const Quadrupole&>(components[compj].quadrupole(qj));
          double absqj=quadrupolej.absQ();
          double q2075=.75*absqi*absqj;
          pstrmij << q2075;
        }
      }
    }
  }
}
