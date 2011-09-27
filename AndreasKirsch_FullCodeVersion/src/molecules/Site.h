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
#include <cmath>

/** Site
 @author Martin Bernreuther et al. (2009)
 */
class Site {
public:
	/// get x-coordinate of position
	double rx() const { return _r[0]; }
	/// get y-coordinate of position
	double ry() const { return _r[1]; }
	/// get z-coordinate of position
	double rz() const { return _r[2]; }
	/// get position vector
	const double* r() const { return _r; }
	/// get mass
	double m() const { return _m; }
	/// translate coordinates to new origin
	void translateOrigin(double neworigin[3]) {
		for (unsigned short d = 0; d < 3; ++d)
			_r[d] -= neworigin[d];
	}

protected:
	/// Constructor
	Site(double x = 0., double y = 0., double z = 0., double m = 0.)
			: _m(m) {
		_r[0] = x;
		_r[1] = y;
		_r[2] = z;
	}

	//private:
	double _r[3]; // position coordinates
	double _m; // mass
};


/** LJcenter
 Lennard-Jones 12-6 center
 @author Martin Bernreuther et al.
 */
class LJcenter : public Site {
public:
	/// Constructor: pass shift = 0.0 for the full LJ potential
	LJcenter(double x, double y, double z, double m, double eps, double sigma, double rc, double shift)
			: Site(x, y, z, m), _eps(eps), _sigma(sigma), _rc(rc) {
		_r[0] = x;
		_r[1] = y;
		_r[2] = z;
		_uLJshift6 = shift;
	}
	/// Constructor reading from stream
	LJcenter(std::istream& istrm) {
		istrm >> _r[0] >> _r[1] >> _r[2] >> _m >> _eps >> _sigma >> _uLJshift6;
	}
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _m << "\t" << _eps << " " << _sigma << " " << _rc << " " << _uLJshift6;
	}
	/// get strength
	double eps() const { return _eps; }
	/// get diameter
	double sigma() const { return _sigma; }
	/// get truncation option
	bool TRUNCATED_SHIFTED() {
		return (this->_uLJshift6 != 0.0);
	}
	double shift6() const {
		return this->_uLJshift6;
	}

private:
	double _eps; // strength
	double _sigma; // diameter
  // cutoff radius
  // it seems to me as if this is not the cutoff-radius which is used for the linked cells,
  // but a molecule specific one to determine the cutoff correction
  // TODO why may they be different!!!???
	double _rc;
	double _uLJshift6; // truncation offset of the LJ potential
};

class Charge : public Site {
public:
	/// Constructor
	Charge(double x, double y, double z, double m, double q)
			: Site(x, y, z, m), _q(q) {
		_r[0] = x; _r[1] = y; _r[2] = z;
		_m = m;
		_q = q;
	}

	/// Constructor reading from stream
	Charge(std::istream& istrm) {
		istrm >> _r[0] >> _r[1] >> _r[2] >> _m >> _q;
	}
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _m << " " << _q;
	}
	/// get charge
	double q() const { return _q; }

private:
	/// charge
	double _q;
};

class Tersoff : public Site {
public:
	Tersoff(double x, double y, double z,
	        double m, double A, double B,
	        double lambda, double mu, double R, double S,
	        double c, double d, double h, double n, double beta)
		: Site(x, y, z, m)
	{
		this->_r[0] = x;
		this->_r[1] = y;
		this->_r[2] = z;

		this->_m = m;
		this->_A = A;
		this->_B = B;

		this->_minus_lambda = -lambda;
		this->_minus_mu = -mu;
		this->_R = R;
		this->_S = S;

		this->_c_square = c*c;
		this->_d_square = d*d;
		this->_h = h;
		this->_n = n;
		this->_beta = beta;
	}

	/// Constructor reading from stream
	Tersoff(std::istream& istrm) {
		double lambda, mu, c, d;
		istrm >> _r[0] >> _r[1] >> _r[2]
		      >> _m >> _A >> _B >> lambda >> mu
		      >> _R >> _S >> c
		      >> d >> _h >> _n >> _beta;
		_minus_lambda = -lambda;
		_minus_mu = -mu;
		_c_square = c*c;
		_d_square = d*d;
	}

	//! @brief write to stream
	//!
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t"
		      << _m << "\t"  << _A << " "  << _B << " " << -_minus_lambda << " " << -_minus_mu << "\t"
		      << _R << " " << _S << "\t" << sqrt(_c_square) << " "
		      << sqrt(_d_square) << "\t" << _h << " " << _n << " " << _beta;
	}

	//! @brief get repulsive coefficient
	//!
	double A() const { return this->_A; }

	//! @brief get attractive coefficient
	//!
	double B() const { return this->_B; }

	//! @brief get repulsive coexponent
	//!
	double minusLambda() const { return this->_minus_lambda; }

	//! @brief get attractive coexponent
	//!
	double minusMu() const { return this->_minus_mu; }

	//! @brief get internal radius
	//!
	double R() const { return this->_R; }

	//! @brief get external radius
	//!
	double S() const { return this->_S; }

	//! @brief get c square Tersoff interaction parameter
	//!
	double cSquare() const { return this->_c_square; }

	//! @brief get d square Tersoff interaction parameter
	//!
	double dSquare() const { return this->_d_square; }

	//! @brief get h Tersoff interaction parameter
	//!
	double h() const { return this->_h; }

	//! @brief get n Tersoff interaction parameter
	//!
	double n() const { return this->_n; }

	//! @brief get beta Tersoff interaction parameter
	//!
	double beta() const { return this->_beta; }

private:
	//! repulsive coefficient
	double _A;
	//! attractive coefficient
	double _B;
	//! repulsive coexponent
	double _minus_lambda;
	//! attractive coexponent
	double _minus_mu;

	//! internal radius
	double _R;
	//! external radius
	double _S;

	//! c square Tersoff interaction parameter
	double _c_square;
	//! d square Tersoff interaction parameter
	double _d_square;
	//! h Tersoff interaction parameter
	double _h;
	//! n Tersoff interaction parameter
	double _n;
	//! beta Tersoff interaction parameter
	double _beta;
};

/** OrientedSite
 @author Martin Bernreuther
 */
class OrientedSite : public Site {
public:
	double ex() const { return _e[0]; }
	double ey() const { return _e[1]; }
	double ez() const { return _e[2]; }
	const double* e() const { return _e; }

protected:
	/// Constructor
	OrientedSite(double x = 0., double y = 0., double z = 0., double m = 0., double ex = 0., double ey = 0., double ez = 0.)
			: Site(x, y, z, m) {
		_e[0] = ex;
		_e[1] = ey;
		_e[2] = ez;
	}

	//private:
	double _e[3];
};


/** Dipole
 @author Martin Bernreuther
 */
class Dipole : public OrientedSite {
public:
	/// Constructor
	Dipole(double x, double y, double z, double eMyx, double eMyy, double eMyz, double absMy)
			: OrientedSite(x, y, z, 0., eMyx, eMyy, eMyz), _absMy(absMy) {
	}
	/// Constructor reading from stream
	Dipole(std::istream& istrm) {
		istrm >> _r[0] >> _r[1] >> _r[2] >> _e[0] >> _e[1] >> _e[2] >> _absMy;
		_m = 0.;
	}
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _e[0] << " " << _e[1] << " " << _e[2] << "\t" << _absMy;
	}
	double absMy() const { return _absMy; }

private:
	double _absMy;
};

/** Quadrupole
 @author Martin Bernreuther
 */
class Quadrupole : public OrientedSite {
public:
	/// Constructor
	Quadrupole(double x, double y, double z, double eQx, double eQy, double eQz, double absQ)
			: OrientedSite(x, y, z, 0., eQx, eQy, eQz), _absQ(absQ) {
	}
	/// Constructor reading from stream
	Quadrupole(std::istream& istrm) {
		istrm >> _r[0] >> _r[1] >> _r[2] >> _e[0] >> _e[1] >> _e[2] >> _absQ;
		_m = 0.;
	}
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _e[0] << " " << _e[1] << " " << _e[2] << " " << _absQ;
	}
	double absQ() const { return _absQ; }

private:
	double _absQ;
};

#endif /*SITE_H_*/
