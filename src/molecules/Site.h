#ifndef SITE_H_
#define SITE_H_

#include <iostream>
#include <cmath>
#include <cassert>

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

using Log::global_log;

/** Site
 *
 * Sites are the center of a interactions. Depending on the interaction type
 * there are specialised derived classes of this basic site class.
 *
 * @author Martin Bernreuther et al. (2009)
 */
class Site {
public:
	double rx() const { return _r[0]; }  /**< get x-coordinate of position vector */
	double ry() const { return _r[1]; }  /**< get y-coordinate of position vector */
	double rz() const { return _r[2]; }  /**< get z-coordinate of position vector */
	const double* r() const { return _r; }  /**< get position vector */
	double m() const { return _m; }         /**< get mass */

	/**
	 * set the d-th component of the position
	 */
	void setR(int d, double r) {
		assert(d < 3);
		_r[d] = r;
	}

	void setM(double m) {
		_m = m;
	}

	virtual ~Site() {}

protected:
	/// Constructor
    Site(double x = 0., double y = 0., double z = 0., double m = 0.)
        : _m(m) {
            _r[0] = x;
            _r[1] = y;
            _r[2] = z;
        }

	double _r[3]; /**< position coordinates */
	double _m;    /**< mass */
};


/** Lennard-Jones 12-6 center
 * @author Martin Bernreuther et al.
 *
 * Lennard-Jones 12-6 interaction site. The potential between two LJ centers of the same type is 
 * given by 
 *
 * \f$[
 *  U_\text{LJ} = \epsilon \left[ \left(\frac{r}{\sigma}\right)^{6} - \left(\frac{r}{\sigma}\right)^{12} \right]
 * \f]
 *
 * where $r$ is the distance between the two LJ centers. See potforce.h for the detailed implementation.
 */
class LJcenter : public Site {
public:
	LJcenter() {}
	/** Constructor
     *
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] m        mass    
     * \param[in] epsilon  interaction strength
     * \param[in] sigma    interaction diameter
     * \param[in] rc       cutoff radius
     * \param[in] shift    0. for full LJ potential
     */
	LJcenter(double x, double y, double z, double m, double eps, double sigma, double rc, double shift)
			: Site(x, y, z, m), _eps(eps), _sigma(sigma), _rc(rc), _uLJshift6(shift) { }

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("mass", _m);
		xmlconfig.getNodeValueReduced("epsilon", _eps);
		xmlconfig.getNodeValueReduced("sigma", _sigma);
		if( xmlconfig.getNodeValueReduced("cutoff", _rc) == 0 ) {
			/* TODO: remove rc from LJ site? */
			_rc = -1.; /* set to invalid value  */
			global_log->warning() << "Cutoff radius for LJ site not specified" << std::endl;
		}
		_uLJshift6 = 0.0;
		xmlconfig.getNodeValueReduced("shifted", _uLJshift6);
	}
	
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _m << "\t" << _eps << " " << _sigma << " " << _rc << " " << _uLJshift6;
	}
    /* TODO rename to epsilon */
	double eps() const { return _eps; }     /**< get interaction strength */
	double sigma() const { return _sigma; } /**< get interaction diameter */

	void setEps(double eps) {
		_eps = eps;
	}

	void setSigma(double sigma) {
		_sigma = sigma;
	}

	void setRC(double rc) {
		_rc = rc;
	}

	void setULJShift6(double uLJshift6) {
		_uLJshift6 = uLJshift6;
	}

    /* TODO: The following method is never used */
	bool TRUNCATED_SHIFTED() { return (_uLJshift6 != 0.0); } /**< get truncation option */

	double shift6() const { return _uLJshift6; }             /**< get energy shift of interaction potential */

private:
	double _eps;    /**< interaction strength */
	double _sigma;  /**< interaction diameter */

	// cutoff radius
	// it seems to me as if this is not the cutoff-radius which is used for the linked cells,
	// but a molecule specific one to determine the cutoff correction
	// TODO why may they be different!!!???
	double _rc;     /**< cutoff radius */
	double _uLJshift6; /**< truncation offset of the LJ potential */
};

/** Charge center
 */
class Charge : public Site {
public:
	Charge() {}
    /** Constructor
     *
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] m        mass    
     * \param[in] q        charge
     */
    Charge(double x, double y, double z, double m, double q)
			: Site(x, y, z, m), _q(q) { }

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("mass", _m);
		xmlconfig.getNodeValueReduced("charge", _q);
	}
	
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _m << " " << _q;
	}
	double q() const { return _q; }  /**< get charge */

	void setQ(double q) {
		_q = q;
	}

private:
	double _q;  /**< charge */
};



/** Tersoff center
 *
 * TODO: Document the potential, in particular the parameters
 */
class Tersoff : public Site {
public:
	Tersoff() {}
	Tersoff(double x, double y, double z,
	        double m, double A, double B,
	        double lambda, double mu, double R, double S,
	        double c, double d, double h, double n, double beta)
		: Site(x, y, z, m), _A(A), _B(B), _minus_lambda(-lambda), _minus_mu(-mu), _R(R), _S(S),
        _c_square(c*c), _d_square(d*d), _h(h),_n(n), _beta(beta) {}


	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("A", _A);
		xmlconfig.getNodeValueReduced("B", _B);
		double lambda, mu;
		xmlconfig.getNodeValueReduced("lambda", lambda); _minus_lambda = -lambda;
		xmlconfig.getNodeValueReduced("mu", mu); _minus_mu = -mu;
		xmlconfig.getNodeValueReduced("R", _R);
		xmlconfig.getNodeValueReduced("S", _S);
		double c, d;
		xmlconfig.getNodeValueReduced("c", c); _c_square = c*c;
		xmlconfig.getNodeValueReduced("d", d); _d_square = d*d;
		xmlconfig.getNodeValueReduced("h", _h);
		xmlconfig.getNodeValueReduced("beta", _beta);
	}
	
	//! @brief write to stream
	//!
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t"
		      << _m << "\t"  << _A << " "  << _B << " " << -_minus_lambda << " " << -_minus_mu << "\t"
		      << _R << " " << _S << "\t" << sqrt(_c_square) << " "
		      << sqrt(_d_square) << "\t" << _h << " " << _n << " " << _beta;
	}

	double A() const { return _A; }  /**< get repulsive coefficient */
	double B() const { return _B; }  /**< get attractive coefficient */
	double minusLambda() const { return _minus_lambda; }  /**< get repulsive coexponent */
	double minusMu() const { return _minus_mu; }  /**< get attractive coexponent */

	double R() const { return _R; }  /**< get internal radius */
	double S() const { return _S; }  /**< get external radius */

	double cSquare() const { return _c_square; }  /**< get c square Tersoff interaction parameter */
	double dSquare() const { return _d_square; }  /**< get d square Tersoff interaction parameter */
	double h() const { return _h; }  /**< get h Tersoff interaction parameter */
	double n() const { return _n; }  /**< get n Tersoff interaction parameter */
	double beta() const { return _beta; }  /**< get beta Tersoff interaction parameter */

private:
	double _A;  /**< repulsive coefficient */
	double _B;  /**< attractive coefficient */
	double _minus_lambda;  /**< repulsive coexponent */
	double _minus_mu;      /**< attractive coexponent */

	double _R;  /**< internal radius */
	double _S;  /**< external radius */

	double _c_square;  /**< c square Tersoff interaction parameter */
	double _d_square;  /**< d square Tersoff interaction parameter */
	double _h;  /**< h Tersoff interaction parameter */
	double _n;  /**< n Tersoff interaction parameter */
	double _beta;  /**< beta Tersoff interaction parameter */
};

/** Oriented site
 * @author Martin Bernreuther
 */
class OrientedSite : public Site {
public:
	double ex() const { return _e[0]; }
	double ey() const { return _e[1]; }
	double ez() const { return _e[2]; }
	const double* e() const { return _e; }  /**< Get pointer to the normalized orientation vector. */

	/** set the d-th component of the orientation vector */
	void setE(int d, double e) {
		assert(d < 3);
		_e[d] = e;
	}

protected:
	/// Constructor
	OrientedSite(double x = 0., double y = 0., double z = 0., double m = 0., double ex = 0., double ey = 0., double ez = 0.)
			: Site(x, y, z, m) {
		_e[0] = ex;
		_e[1] = ey;
		_e[2] = ez;
	}

	double _e[3];  /**< Normalized orientation vector */
};


/** Dipole
 * @author Martin Bernreuther
 *
 */
class Dipole : public OrientedSite {
public:
	Dipole() {}
    /** Constructor
     *
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] eMyx     x coordinate of the dipole moments normal
     * \param[in] eMyy     y coordinate of the dipole moments normal
     * \param[in] eMyz     z coordinate of the dipole moments normal
     * \param[in] absQ     dipole moments absolute value
     */
    Dipole(double x, double y, double z, double eMyx, double eMyy, double eMyz, double absMy)
			: OrientedSite(x, y, z, 0., eMyx, eMyy, eMyz), _absMy(absMy) {
	}

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("dipolemoment/abs", _absMy);
		xmlconfig.getNodeValueReduced("dipolemoment/x", _e[0]);
		xmlconfig.getNodeValueReduced("dipolemoment/y", _e[1]);
		xmlconfig.getNodeValueReduced("dipolemoment/z", _e[2]);
		/* TODO normalization check */
	}

	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _e[0] << " " << _e[1] << " " << _e[2] << "\t" << _absMy;
	}

	double absMy() const { return _absMy; }  /**< Get the absolute value of the dipole moment. */

	/** set the value of the dipole moment */
	void setAbyMy(double my) {
		_absMy = my;
	}

private:
    /* TODO: move abs to oriented site. */
	double _absMy;  /**< absolute value of the dipole moment. */
};

/** Quadrupole
 * @author Martin Bernreuther
 */
class Quadrupole : public OrientedSite {
public:
	Quadrupole() {}
    /** Constructor
     *
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] eQx      x coordinate of the quadrupole moments normal
     * \param[in] eQy      y coordinate of the quadrupole moments normal
     * \param[in] eQz      z coordinate of the quadrupole moments normal
     * \param[in] absQ     quadrupole moments absolute value
     */
	Quadrupole(double x, double y, double z, double eQx, double eQy, double eQz, double absQ)
			: OrientedSite(x, y, z, 0., eQx, eQy, eQz), _absQ(absQ) {
	}

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("quadrupolemoment/abs", _absQ);
		xmlconfig.getNodeValueReduced("quadrupolemoment/x", _e[0]);
		xmlconfig.getNodeValueReduced("quadrupolemoment/y", _e[1]);
		xmlconfig.getNodeValueReduced("quadrupolemoment/z", _e[2]);
		/* TODO normalization check */
	}
	
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _e[0] << " " << _e[1] << " " << _e[2] << " " << _absQ;
	}
	double absQ() const { return _absQ; }  /**< Get the absolute value of the quadrupole moment. */

	/** set the absolute value of teh quadrupole moment */
	void setAbsQ(double q) {
		_absQ = q;
	}

private:
    /* TODO: move abs to oriented site. */
	double _absQ;  /**< absolute value of the quadrupole moment. */
};

#endif  /* SITE_H_ */
