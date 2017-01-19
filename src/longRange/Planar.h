/***************************************************************************
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
 *                                                                         *
 *   Due to copyleft all future versions of this program must be           *
 *   distributed as Free Software (e.g., using a BSD-like license).        *
 ***************************************************************************/

//Calculation of the surface tension in a system with planar interfaces needs a Long Range Correction.
//
//The correction terms are based on Janecek (2006) and Lustig (1988).

#ifndef PLANAR_H_
#define PLANAR_H_

#include "LongRangeCorrection.h"

#include <vector>
#include <cmath>
#include <string>
#include <map>

class Domain;
class ParticleContainer;
class Molecule;

class Planar:public LongRangeCorrection {
public:
	Planar(double cutoffT,double cutoffLJ,Domain* domain,  DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, unsigned slabs, Simulation _simulation);
//	 ~Planar();	

	void calculateLongRange();
	double lrcLJ(Molecule* mol);
	void SetSmoothDensityProfileOption(bool bVal) {_smooth = bVal;}

private:

	void centerCenter(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj,double centerM); 
	void centerSite(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj);
	void siteSite(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj);
	void dipoleDipole(unsigned ci,unsigned cj,unsigned si,unsigned sj);

	unsigned _slabs;
	unsigned numComp;
	unsigned *numLJ;
	unsigned *numDipole;
	unsigned *numCharge;
	unsigned numLJSum;
	unsigned numDipoleSum;
	unsigned *numLJSum2;
	unsigned *numDipoleSum2;
	bool _smooth;
	bool _dipole;
	double *uLJ;
	double *vNLJ;
	double *vTLJ;
	double *fLJ;
	double *rho_g;
	double *rho_l;
	double *fDipole;
	double *uDipole;
	double *vNDipole;
	double *vTDipole;
	double *rhoDipole;
	double *rhoDipoleL;
	double *muSquare;
	double *eLong;
	double cutoff;
	double delta;
	unsigned cutoff_slabs;
	int frequency;
	double ymax;
	double boxlength[3];
	double V;
	int sint;
	double temp;
	unsigned simstep;
	
	ParticleContainer* _particleContainer;
	Domain* _domain;
	DomainDecompBase* _domainDecomposition;
	
};


#endif /*Planar_H_*/
