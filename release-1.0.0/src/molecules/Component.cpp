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
#include <iostream>
#include <iomanip>

#include "Component.h"
#include "Site.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

Component::Component(unsigned int id) {
	_id = id;
	_m = 0.;
	_I[0] = _I[1] = _I[2] = _I[3] = _I[4] = _I[5] = 0.;
	_rot_dof = 0;
	_Ipa[0] = _Ipa[1] = _Ipa[2] = 0.;
	_numMolecules = 0;
	this->maximalTersoffExternalRadius = 0.0;

	_ljcenters = vector<LJcenter> ();
	_charges = vector<Charge> ();
	_quadrupoles = vector<Quadrupole> ();
	_dipoles = vector<Dipole> ();
	_tersoff = vector<Tersoff> ();
}

void Component::readXML(XMLfileUnits& xmlconfig) {
	global_log->info() << "Reading in component" << endl;
	unsigned int cid = 0;
	xmlconfig.getNodeValue( "@id", cid );
	global_log->info() << "Component ID:" << cid << endl;
	setID(cid - 1);
	string name;
	xmlconfig.getNodeValue( "@name", name );
	global_log->info() << "Component name:" << name << endl;
	setName(name);

	XMLfile::Query query = xmlconfig.query( "site" );
	XMLfile::Query::const_iterator siteIter;
	for( siteIter = query.begin(); siteIter; siteIter++ ) {
		xmlconfig.changecurrentnode( siteIter );
		
		std::string siteType;
		xmlconfig.getNodeValue( "@type", siteType );
		std::cout << "Adding site of type " << siteType << endl;
		
		if ( siteType == "LJ126" ) {
			LJcenter ljSite;
			ljSite.readXML(xmlconfig);
			addLJcenter( ljSite );
		} else
		if ( siteType == "Charge" ) {
			Charge chargeSite;
			chargeSite.readXML(xmlconfig);
			addCharge(chargeSite);
		} else
		if ( siteType == "Dipole" ) {
			Dipole dipoleSite;
			dipoleSite.readXML(xmlconfig);
			addDipole(dipoleSite);
		} else
		if ( siteType == "Quadrupole" ) {
			Quadrupole quadrupoleSite;
			quadrupoleSite.readXML(xmlconfig);
			addQuadrupole(quadrupoleSite);
		} else
		if ( siteType == "Tersoff" ) {
			Tersoff tersoffSite;
			tersoffSite.readXML(xmlconfig);
			addTersoff(tersoffSite);
		}
	}

	if(xmlconfig.changecurrentnode("momentsofinertia")){
		double II[3];
		if(xmlconfig.getNodeValueReduced("Ixx", II[0]) > 0) { setI11(II[0]); }
		if(xmlconfig.getNodeValueReduced("Iyy", II[1]) > 0) { setI22(II[1]); }
		if(xmlconfig.getNodeValueReduced("Izz", II[2]) > 0) { setI33(II[2]); }
		xmlconfig.changecurrentnode("..");
	}
}



void Component::addLJcenter(double x, double y, double z,
                            double m, double eps, double sigma,
                            double rc, bool TRUNCATED_SHIFTED) {
	double shift6 = 0.0;
	if (TRUNCATED_SHIFTED) {
		double sigperrc2 = sigma * sigma / (rc * rc);
		double sigperrc6 = sigperrc2 * sigperrc2 * sigperrc2;
		shift6 = 24.0 * eps * (sigperrc6 - sigperrc6 * sigperrc6);
	}

	LJcenter ljsite(x, y, z, m, eps, sigma, rc, shift6);
	_ljcenters.push_back(ljsite);
	updateMassInertia(ljsite);
}

void Component::addLJcenter(LJcenter& ljsite) {
	_ljcenters.push_back(ljsite);
	updateMassInertia(ljsite);
}


void Component::updateMassInertia() {
	_m = 0;
	for (int i = 0; i < 6; i++) {
		_I[i] = 0.0;
	}

	for (size_t i = 0; i < _ljcenters.size(); i++) {
		updateMassInertia(_ljcenters[i]);
	}
	for (size_t i = 0; i < _charges.size(); i++) {
		updateMassInertia(_charges[i]);
	}
	for (size_t i = 0; i < _tersoff.size(); i++) {
		updateMassInertia(_tersoff[i]);
	}
}

void Component::updateMassInertia(Site& site) {
	_m += site.m();
	// assume the input is already transformed to the principal axes system
	// (and therefore the origin is the center of mass)
//	_I[0] += m * (y * y + z * z);
	_I[0] += site.m() * (site.ry() * site.ry() + site.rz() * site.rz());
//	_I[1] += m * (x * x + z * z);
	_I[1] += site.m() * (site.rx() * site.rx() + site.rz() * site.rz());
//	_I[2] += m * (x * x + y * y);
	_I[2] += site.m() * (site.rx() * site.rx() + site.ry() * site.ry());
//	_I[3] -= m * x * y;
	_I[3] -= site.m() * site.rx() * site.ry();
//	_I[4] -= m * x * z;
	_I[4] -= site.m() * site.rx() * site.rz();
//	_I[5] -= m * y * z;
	_I[5] -= site.m() * site.ry() * site.rz();

	_rot_dof = 3;
	for (unsigned short d = 0; d < 3; ++d) {
		_Ipa[d] = _I[d];
		if (_Ipa[d] == 0.) --_rot_dof;
	}
}

void Component::addCharge(double x, double y, double z, double m, double q) {
	Charge chargesite(x, y, z, m, q);
	_charges.push_back(chargesite);
	updateMassInertia(chargesite);
}

void Component::addCharge(Charge& chargesite) {
	_charges.push_back(chargesite);
	updateMassInertia(chargesite);
}


void Component::addDipole(double x, double y, double z,
                          double eMyx, double eMyy, double eMyz, double eMyabs) {
	_dipoles.push_back(Dipole(x, y, z, eMyx, eMyy, eMyz, eMyabs));
	// massless...
}

void Component::addDipole(Dipole& dipolesite) {
	_dipoles.push_back(dipolesite);
	// massless...
}


void Component::addQuadrupole(double x, double y, double z,
                              double eQx, double eQy, double eQz, double eQabs) {
	_quadrupoles.push_back(Quadrupole(x, y, z, eQx, eQy, eQz, eQabs));
	// massless...
}

void Component::addQuadrupole(Quadrupole& quadrupolesite) {
	_quadrupoles.push_back(quadrupolesite);
	// massless
}


void Component::addTersoff(double x, double y, double z,
                           double m, double A, double B, double lambda, double mu, double R,
                           double S, double c, double d, double h, double n, double beta) {
	if (S > this->maximalTersoffExternalRadius) maximalTersoffExternalRadius = S;
	Tersoff tersoffsite(x, y, z, m, A, B, lambda, mu, R, S, c, d, h, n, beta);
	_tersoff.push_back(tersoffsite);
	updateMassInertia(tersoffsite);
}

void Component::addTersoff(Tersoff& tersoffsite)
{
	_tersoff.push_back(tersoffsite);
	updateMassInertia(tersoffsite);
}


void Component::write(std::ostream& ostrm) const {
	ostrm << _ljcenters.size() << "\t" << _charges.size() << "\t"
	      << _dipoles.size() << "\t" << _quadrupoles.size() << "\t"
	      << _tersoff.size() << "\n";
	for (std::vector<LJcenter>::const_iterator pos = _ljcenters.begin(); pos != _ljcenters.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	for (std::vector<Charge>::const_iterator pos = _charges.begin(); pos != _charges.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	for (std::vector<Dipole>::const_iterator pos = _dipoles.begin(); pos != _dipoles.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	for (std::vector<Quadrupole>::const_iterator pos = _quadrupoles.begin(); pos != _quadrupoles.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	for (std::vector<Tersoff>::const_iterator pos = _tersoff.begin(); pos != _tersoff.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	ostrm << _Ipa[0] << " " << _Ipa[1] << " " << _Ipa[2] << endl;
}

void Component::writePOVobjs(std::ostream& ostrm, string para) const {
	if (numLJcenters() <= 0) return;
	if (numLJcenters() == 1) {
		ostrm << "sphere {<" << _ljcenters.front().rx() << "," << _ljcenters.front().ry() << "," << _ljcenters.front().rz() << ">," << .5 * _ljcenters.front().sigma() << " " << para << "}";
	}
	else {
		ostrm << "blob { threshold 0.01 ";
		for (std::vector<LJcenter>::const_iterator pos = _ljcenters.begin(); pos != _ljcenters.end(); ++pos)
			ostrm << "sphere {<" << pos->rx() << "," << pos->ry() << "," << pos->rz() << ">," << .5 * pos->sigma() << ", strength 1 } ";
		ostrm << para << "}";
	}
	ostrm << flush;
}

void Component::writeVIM(std::ostream& ostrm) {
	for (std::vector<LJcenter>::const_iterator pos = _ljcenters.begin(); pos != _ljcenters.end(); ++pos) {
		ostrm << "~ " << this->_id + 1 << " LJ " << setw(7) << pos->rx() << ' '
		      << setw(7) << pos->ry() << ' ' << setw(7) << pos->rz() << ' '
		      << setw(6) << pos->sigma() << ' ' << setw(2) << (1 + (this->_id % 9)) << "\n";
	}
	ostrm << flush;
}

std::ostream& operator<<(std::ostream& stream, const Component& component) {
	component.write(stream);
	return stream;
}
