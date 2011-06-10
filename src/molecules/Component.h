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

#ifndef COMPONENT_H_
#define COMPONENT_H_

#include <vector>
#include <string>

#include "molecules/Site.h"

/**
 Komponente (Molekueltyp)

 @author Martin Bernreuther
 */
class Component {
public:
	Component(unsigned int id = 0);

	void setID(unsigned int id) { _id = id; }
	unsigned int ID() const { return _id; }
	unsigned int numSites() const {
		return this->numLJcenters() + this->numCharges()
		                            + this->numDipoles()
		                            + this->numQuadrupoles()
		                            + this->numTersoff();
	}
	unsigned int numLJcenters() const { return _ljcenters.size(); }
	unsigned int numCharges() const { return _charges.size(); }
	unsigned int numDipoles() const { return _dipoles.size(); }
	unsigned int numQuadrupoles() const { return _quadrupoles.size(); }
	unsigned int numTersoff() const { return _tersoff.size(); }

	double m() const { return _m; }
	double I11() const { return _Ipa[0]; }
	double I22() const { return _Ipa[1]; }
	double I33() const { return _Ipa[2]; }
	void setI11(double I) { _Ipa[0]=I; }
	void setI22(double I) { _Ipa[1]=I; }
	void setI33(double I) { _Ipa[2]=I; }

	unsigned int getRotationalDegreesOfFreedom() const { return _rot_dof; }

	const std::vector<LJcenter>& ljcenters() const { return _ljcenters; }
	LJcenter& ljcenter(unsigned int i) { return _ljcenters[i]; }
	const LJcenter& ljcenter(unsigned int i) const { return _ljcenters[i]; }
	const std::vector<Charge>& charges() const { return _charges; }
	Charge& charge(unsigned i) { return _charges[i]; }
	const Charge& charge(unsigned i) const { return _charges[i]; }
	const std::vector<Dipole>& dipoles() const { return _dipoles; }
	Dipole& dipole(unsigned int i) { return _dipoles[i]; }
	const Dipole& dipole(unsigned int i) const { return _dipoles[i]; }
	const std::vector<Quadrupole>& quadrupoles() const { return _quadrupoles; }
	Quadrupole& quadrupole(unsigned int i) { return _quadrupoles[i]; }
	const Quadrupole& quadrupole(unsigned int i) const { return _quadrupoles[i]; }
	const std::vector<Tersoff>& tersoff() const { return _tersoff; }
	Tersoff& tersoff(unsigned int i) { return _tersoff[i]; }
	const Tersoff& tersoff(unsigned int i) const { return _tersoff[i]; }

	void setNumMolecules(unsigned long num) { _numMolecules = num; }  /**< set the number of molecules for this component */
	void incNumMolecules() { ++_numMolecules; }  /**< increase the number of molecules for this component by 1 */ 
	void incNumMolecules(int N) { _numMolecules += N; }  /**< increase the number of molecules for this component by N */
	unsigned long getNumMolecules() const { return _numMolecules; }  /**< get the number of molecules (global) of the component */

	void addLJcenter(
			double x, double y, double z, double m, double eps,
			double sigma, double rc = 0, bool TRUNCATED_SHIFTED = 0
	);
	void addCharge(double x, double y, double z, double m, double q);
	void addDipole(double x, double y, double z,
	               double eMyx, double eMyy, double eMyz, double eMyabs);
	void addQuadrupole(double x, double y, double z,
	                   double eQx, double eQy, double eQz, double eQabs);
	void addTersoff(double x, double y, double z,
	                double m, double A, double B, double lambda, double mu, double R,
	                double S, double c, double d, double h, double n, double beta);

	/** delete the last site stored in the vector */
	void deleteLJCenter() { _ljcenters.pop_back() ;}
	void deleteCharge() { _charges.pop_back() ;}
	void deleteDipole() { _dipoles.pop_back() ;}
	void deleteQuadrupole() { _quadrupoles.pop_back() ;}
	void deleteTersoff() { _tersoff.pop_back() ;}

	/**
	 * To be called after sites have been deleted or the properties of sites have been changed.
	 */
	void updateMassInertia();

	/** write information to stream */
	void write(std::ostream& ostrm) const;

	/** write POVray object definition to stream */
	void writePOVobjs(std::ostream& ostrm, std::string para = std::string("pigment {color rgb<1,0,0>}")) const;

	void writeVIM(std::ostream& ostrm);

	double getTersoffRadius() { return this->maximalTersoffExternalRadius; }
	void setTersoffRadius(double mTER) { this->maximalTersoffExternalRadius = mTER; }

	void setE_trans(double E) { _E_trans = E; }
	void setE_rot(double E) { _E_rot = E; }
	void setT(double T) { _T = T; }
	double E_trans() { return _E_trans; }
	double E_rot() { return _E_rot; }
	double E() { return _E_trans + _E_rot; }
	double T() { return _T; }

private:

	void updateMassInertia(Site& site);

	unsigned int _id; // IDentification number
	// LJcenter,Dipole,Quadrupole have different size -> not suitable to store in a _Site_-array
	//std::vector<Site> _sites;
	// use separate vectors instead...
	std::vector<LJcenter> _ljcenters;
	std::vector<Charge> _charges;
	std::vector<Dipole> _dipoles;
	std::vector<Quadrupole> _quadrupoles;
	std::vector<Tersoff> _tersoff;
	// for performance reasons better(?) omit Site-class indirection and use cached values
	double _m; // total mass
	// Ixx,Iyy,Izz,Ixy,Ixz,Iyz
	double _I[6]; // moments of inertia tensor
	unsigned long _rot_dof; // number of rotational degrees of freedom
	double _Ipa[3]; // moments of inertia for principal axes

	/* cached values, set by ensemble class! */
	unsigned long _numMolecules; // number of molecules for this molecule type
	double _E_trans; // translational energy
	double _E_rot; // rotational energy
	double _T; // temperature

	double maximalTersoffExternalRadius;
};

std::ostream& operator<<(std::ostream& stream, const Component& component);

#endif /*COMPONENT_H_*/
