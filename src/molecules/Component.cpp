#include <iostream>
#include <iomanip>

#include "Component.h"
#include "Site.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "Simulation.h"

using namespace std;
using Log::global_log;

Component::Component(unsigned int id) {
	_id = id;
	_m = 0.;
	_I[0] = _I[1] = _I[2] = _I[3] = _I[4] = _I[5] = 0.;
	_rot_dof = 0;
	_Ipa[0] = _Ipa[1] = _Ipa[2] = 0.;
	_numMolecules = 0;
	_T = 1.;
	_E_trans=0.;
	_lookUpID = 0;
	_E_rot=0.;
	_isStockmayer = false;

	_ljcenters = vector<LJcenter> ();
	_charges = vector<Charge> ();
	_quadrupoles = vector<Quadrupole> ();
	_dipoles = vector<Dipole> ();
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
		xmlconfig.changecurrentnode(siteIter);

		std::string siteType;
		xmlconfig.getNodeValue("@type", siteType);
		global_log->info() << "Adding site of type " << siteType << endl;

		if (siteType == "LJ126") {
			LJcenter ljSite;
			ljSite.readXML(xmlconfig);
			addLJcenter(ljSite);
		} else if (siteType == "Charge") {
			Charge chargeSite;
			chargeSite.readXML(xmlconfig);
			addCharge(chargeSite);
		} else if (siteType == "Dipole") {
			Dipole dipoleSite;
			dipoleSite.readXML(xmlconfig);
			addDipole(dipoleSite);
		} else if (siteType == "Stockmayer") {
			_isStockmayer = true;
			_rot_dof = 2;

			_Ipa[0] = 1.0;
			_Ipa[1] = 1.0;
			_Ipa[2] = 0.0;

			global_log->info() << "Rotation enabled with [Ixx Iyy Izz] = [" << _Ipa[0] << " " << _Ipa[1] << " "
							   << _Ipa[2] << "]. Dipole direction vector of the Stockmayer fluid should be [0 0 1]."
							   << endl;

		} else if (siteType == "Quadrupole") {
			Quadrupole quadrupoleSite;
			quadrupoleSite.readXML(xmlconfig);
			addQuadrupole(quadrupoleSite);
		} else if (siteType == "Tersoff") {
			global_log->error() << "Tersoff no longer supported:" << siteType << endl;
			Simulation::exit(-1);
		} else {
			global_log->error() << "Unknown site type:" << siteType << endl;
			Simulation::exit(-1);
		}
		// go back to initial level, to be consistent, even if no site information is found.
		xmlconfig.changecurrentnode("..");
	}

	if(xmlconfig.changecurrentnode("momentsofinertia")){
		double II[3];
		if(xmlconfig.getNodeValueReduced("Ixx", II[0]) > 0) { setI11(II[0]); }
		if(xmlconfig.getNodeValueReduced("Iyy", II[1]) > 0) { setI22(II[1]); }
		if(xmlconfig.getNodeValueReduced("Izz", II[2]) > 0) { setI33(II[2]); }
		global_log->info() << "Using moments of inertia set in xml config: Ixx = " << I11() << " ; Iyy = " << I22() << " ; Izz = " << I33() << std::endl;
		xmlconfig.changecurrentnode("..");
	} else {
		global_log->info() << "Using calculated moments of inertia: Ixx = " << I11() << " ; Iyy = " << I22() << " ; Izz = " << I33() << std::endl;
	}
}



void Component::addLJcenter(double x, double y, double z,
                            double m, double eps, double sigma,
                            double rc, bool TRUNCATED_SHIFTED) {
	double shift6 = 0.0;
	if (TRUNCATED_SHIFTED) {
		shift6 = calculateLJshift(eps, sigma, rc);
	}

	LJcenter ljsite(x, y, z, m, eps, sigma, shift6);
	_ljcenters.push_back(ljsite);
	updateMassInertia(ljsite);
}

double Component::calculateLJshift(double eps, double sigma, double rc) const {
	const double sigperrc2 = sigma * sigma / (rc * rc);
	const double sigperrc6 = sigperrc2 * sigperrc2 * sigperrc2;
	return 24.0 * eps * (sigperrc6 - sigperrc6 * sigperrc6);
}

void Component::addLJcenter(LJcenter& ljsite) {
	_ljcenters.push_back(ljsite);
	updateMassInertia(ljsite);
}

void Component::updateAllLJcenters(double rc) {
	for(LJcenter &ljcenter : _ljcenters) {
		if(ljcenter.shiftRequested())
			ljcenter.setULJShift6(calculateLJshift(ljcenter.eps(), ljcenter.sigma(), rc));
	}
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
}

void Component::updateMassInertia(Site& site) {
	_m += site.m();
	// assume the input is already transformed to the principal axes system
	// (and therefore the origin is the center of mass)
	
	if ( not _isStockmayer){ //if the component is a Stockmayer fluid, the moments of inertia are fixed at [1 1 0]
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




void Component::write(std::ostream& ostrm) const {
	ostrm << _ljcenters.size() << "\t" << _charges.size() << "\t"
	      << _dipoles.size() << "\t" << _quadrupoles.size() << "\t"
		  << 0 << "\n";  // the 0 indicates a zero amount of tersoff sites.
	for (auto pos = _ljcenters.cbegin(); pos != _ljcenters.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	for (auto pos = _charges.cbegin(); pos != _charges.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	for (auto pos = _dipoles.cbegin(); pos != _dipoles.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	for (auto pos = _quadrupoles.cbegin(); pos != _quadrupoles.end(); ++pos) {
		pos->write(ostrm);
		ostrm << endl;
	}
	ostrm << _Ipa[0] << " " << _Ipa[1] << " " << _Ipa[2] << endl;
}

void Component::writeVIM(std::ostream& ostrm) {
	for (auto pos = _ljcenters.cbegin(); pos != _ljcenters.end(); ++pos) {
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
