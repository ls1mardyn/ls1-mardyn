/*************************************************************************
 * Copyright (C) 2012 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
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

#ifndef HOMOGENEOUS_H__
#define HOMOGENEOUS_H__

#include "LongRangeCorrection.h"

#include <cmath>

class Domain;
class LongRangeCorrection;

//class Homogeneous:public LongRangeCorrection{
class Homogeneous: public LongRangeCorrection{

public:
//	Homogeneous();
	Homogeneous(double cutoffRadius, double cutoffRadiusLJ,  Domain* domain, Simulation _simulation);
  
//	void initializeLongRange();
	void calculateLongRange();

private:
	double _UpotCorr;
	double _VirialCorr;
	/* TODO: Comments on all the functions */
	// Long range correction for the Lennard-Jones interactions based on Lustig (1988)
	double _TICCu(int n,double rc,double sigma2);
	double _TICSu(int n,double rc,double sigma2,double tau);
	double _TISSu(int n,double rc,double sigma2,double tau1,double tau2);
	double _TICCv(int n,double rc,double sigma2);
	double _TICSv(int n,double rc,double sigma2,double tau);
	double _TISSv(int n,double rc,double sigma2,double tau1,double tau2);
	
	//! Components resp. molecule types
	std::vector<Component> _components;
	//! parameter streams for each possible pair of molecule-types
	Comp2Param _comp2params;
	
	Domain* _domain;
};

#endif /* __HOMOGENEOUS_H__ */
