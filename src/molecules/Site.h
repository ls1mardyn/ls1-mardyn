#ifndef SITE_H_
#define SITE_H_

#include "utils/mardyn_assert.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"
#include <array>
#include <cmath>
#include <cstdint>



/** @brief The Site class is the basis for the implementation of physical interactions.
 *
 * A site in ls1-MarDyn is the center of an interaction. Depending on the interaction type
 * there are specialised derived classes of this basic site class.
 * All sites provide as common propery the site position and site mass.
 */
class Site {
public:
	double rx() const { return _r[0]; }  /**< get x-coordinate of position vector */
	double ry() const { return _r[1]; }  /**< get y-coordinate of position vector */
	double rz() const { return _r[2]; }  /**< get z-coordinate of position vector */
	std::array<double, 3> r() const { return _r; }  /**< get position vector */
	double m() const { return _m; }         /**< get mass */

	/** set the d-th component of the position */
	void setR(int d, double r) {
		mardyn_assert(d < 3);
		_r[d] = r;
	}

	/** set the site mass */
	void setM(double m) { _m = m; }

	/** set the site name */
	void setName(std::string& n) { _name = n; }

	/** get the site name */
	std::string getName() const { return _name; }

	/** @brief Read in XML configuration for a site and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <site type="LJ126|Charge|Dipole|Quadrupole" id="UINT" name="STRING">
	     <coords> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </coords>
	     <mass>DOUBLE</mass>
	   </site>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig) {
		if (!xmlconfig.getNodeValue("@name", _name)) {
			Log::global_log->warning() << "Cannot find site name. Defaulting to type." << std::endl;
			xmlconfig.getNodeValue("@type", _name);
		}
		Log::global_log->info() << "Site name: " << _name << std::endl;
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		Log::global_log->info() << "Site coordinates: (x, y, z) = (" << _r[0] << ", " << _r[1] << ", " << _r[2] << ")" << std::endl;
		xmlconfig.getNodeValueReduced("mass", _m);
		Log::global_log->info() << "Site mass: " << _m << std::endl;
	}

	virtual ~Site() {}

	virtual void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2];
	}

protected:

	/** @brief Constructor
	 * @param x  elative x-coordinate
	 * @param y  elative y-coordinate
	 * @param z  elative z-coordinate
	 * @param m  mass
	 */
    Site(double x = 0., double y = 0., double z = 0., double m = 0.) : _r{{x, y, z}}, _m(m)  {}

	std::array<double, 3> _r;  /**< position coordinates */
	double _m;  /**< mass */
	std::string _name;
};


/** @brief Lennard-Jones 12-6 center
 *
 * Lennard-Jones 12-6 interaction site. The potential between two LJ centers of the same type is
 * given by
 * \f[
 *  U_\text{LJ} = \epsilon \left[ \left(\frac{r}{\sigma}\right)^{6} - \left(\frac{r}{\sigma}\right)^{12} \right]
 * \f]
 * where \f$r\f$ is the distance between the two LJ centers. See potforce.h for the detailed implementation.
 */
class LJcenter : public Site {
public:
	/** @brief Constructor */
	LJcenter(): Site(0., 0., 0., 0.), _epsilon(0.), _sigma(0.), _uLJshift6(0.), _shiftRequested(false) {}
	/** @brief Constructor
	 * \param[in] x        relative x coordinate
	 * \param[in] y        relative y coordinate
	 * \param[in] z        relative z coordinate
	 * \param[in] m        mass
	 * \param[in] epsilon  interaction strength
	 * \param[in] sigma    interaction diameter
	 * \param[in] shift    0. for full LJ potential
	 */
	LJcenter(double x, double y, double z, double m, double epsilon, double sigma, double shift)
		: Site(x, y, z, m), _epsilon(epsilon), _sigma(sigma), _uLJshift6(shift), _shiftRequested(false) {}

	/** @brief Read in XML configuration for a LJcenter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <site>
	     <!-- all Site class parameters -->
	     <epsilon>DOUBLE</epsilon>
	     <sigma>DOUBLE</sigma>
	     <shifted>BOOLEAN</shifted>
	   </site>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) {
		Site::readXML(xmlconfig);
		Log::global_log->info() << "Site type: Lennard-Jones 12-6" << std::endl;
		xmlconfig.getNodeValueReduced("epsilon", _epsilon);
		xmlconfig.getNodeValueReduced("sigma", _sigma);
		xmlconfig.getNodeValue("shifted", _shiftRequested);
		Log::global_log->info() << "Site parameters: epsilon = " << _epsilon << ", sigma: " << _sigma << ", shifted: " << _shiftRequested << std::endl;
	}

	/// write to stream
	void write(std::ostream& ostrm) const {
		Site::write(ostrm);
		ostrm << "\t" << eps() << " " << sigma() << " " << shift6();
	}
    /* @todo TODO: rename to epsilon */
	double eps() const { return _epsilon; }  /**< get interaction strength */
	double sigma() const { return _sigma; }  /**< get interaction diameter */
	double shift6() const { return _uLJshift6; }  /**< get energy shift of interaction potential */
	bool shiftRequested() const { return _shiftRequested; } /**< get the shift request value */

	/** set the interaction strength */
	void setEps(double epsilon) { _epsilon = epsilon; }
	/** set the interaction diameter */
	void setSigma(double sigma) { _sigma = sigma; }
	/** set the energy shift of the interaction potential */
	void setULJShift6(double uLJshift6) { _uLJshift6 = uLJshift6; }


private:
	double _epsilon;  /**< interaction strength */
	double _sigma;  /**< interaction diameter */
	double _uLJshift6; /**< energy shift of the interaction potential, used to implement the LJ truncated and shifted (LJTS) potential */
	bool _shiftRequested; /***< whether the LJTS potential shift needs to be calculated or not */
};

/** @brief Charge center
 *
 * Electrical charge interaction site. The potential between two charge centers is
 * given by the coulomb interaction
 * \f[
 *  U_\text{coulomb} = \frac{1}{4 \pi \epsilon_0} \frac{q_1 q_2}{r}
 * \f]
 * where \f$r\f$ is the distance between the two coulomb centers, \f$q_1\f$ and
 * \f$q_2\f$ are the charges of the two centers and \f$\epsilon_0\f$ is the vacuum
 * permittivity. See potforce.h for the detailed implementation.
 */
class Charge : public Site {
public:
	Charge() : Site(0., 0., 0., 0.), _q(0.) {}
	/** @brief Constructor
	 * \param[in] x        relative x coordinate
	 * \param[in] y        relative y coordinate
	 * \param[in] z        relative z coordinate
	 * \param[in] m        mass
	 * \param[in] q        charge
	 */
    Charge(double x, double y, double z, double m, double q) : Site(x, y, z, m), _q(q) {}


	/** @brief Read in XML configuration for a Charge and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <site>
	     <!-- all Site class parameters -->
	     <charge>DOUBLE</charge>
	   </site>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) {
		Site::readXML(xmlconfig);
		Log::global_log->info() << "Site type: Charge" << std::endl;
		xmlconfig.getNodeValueReduced("charge", _q);
		Log::global_log->info() << "Site parameters: charge = " << _q << std::endl;
	}

	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _m << " " << _q;
	}

	double q() const { return _q; }  /**< get charge */

	/** set the charge */
	void setQ(double q) { _q = q; }

private:
	double _q;  /**< charge */
};


/** @brief The OrientedSite class is the basis for the implementation of physical interactions
 *         based on directed quantities like the dipole moment, etc.
 */
class OrientedSite : public Site {
public:
	double ex() const { return _e[0]; }  /**< get x-coordinate of the normalized orientation vector */
	double ey() const { return _e[1]; }  /**< get y-coordinate of the normalized orientation vector */
	double ez() const { return _e[2]; }  /**< get z-coordinate of the normalized orientation vector */
	std::array<double, 3> e() const { return _e; }  /**< get normalized orientation vector. */
	double abs() const { return _abs; }  /**< get the absolute value of the directed quantity */


	/// write to stream
	void write(std::ostream& ostrm) const {
		Site::write(ostrm);
		ostrm << "\t" << ex() << " " << ey() << " " << ez() << "\t" << abs();
	}

	/** @brief set orientation vector using polar angles
	 * @param theta_deg  theta in degrees
	 * @param phi_deg    phi in degrees
	 */
	void setOrientationVectorByPolarAngles(const double& theta_deg, const double& phi_deg) {
		const double fac = M_PI / 180.;  // translate: degrees --> rad
		double theta_rad = theta_deg * fac;
		double phi_rad   = phi_deg   * fac;
		_e[0] = sin(theta_rad) * cos(phi_rad);
		_e[1] = sin(theta_rad) * sin(phi_rad);
		_e[2] = cos(theta_rad);
	}
	/** set the absolute value of the directed quantity */
	void setAbs(double abs) {_abs = abs; }

	/** set the d-th component of the orientation vector
	 * @note after setting components of the orientation vector normalization has to be performed
	 *       manually using the normalize_e() method.
	 */
	void setE(int d, double e) {
		mardyn_assert(d < 3);
		_e[d] = e;
	}
	/** normalize the orientation vector */
	void normalize_e() {
		_abs = sqrt(ex()*ex() + ey()*ey() + ez()*ez());
		for(int d = 0; d < 3; ++d) { _e[d] /= _abs; }
	}

protected:
    /** @brief Constructor
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] ex       x coordinate of the orientation normal
     * \param[in] ey       y coordinate of the orientation normal
     * \param[in] ez       z coordinate of the orientation normal
     * \param[in] abs      absolute value of the directed quantity
     */
	OrientedSite(double x = 0., double y = 0., double z = 0., double m = 0., double ex = 0., double ey = 0., double ez = 0., double abs = 0.)
		: Site(x, y, z, m), _e{{ex, ey, ez}}, _abs(abs) {}

	std::array<double, 3> _e;  /**< orientation vector @todo normalized!? */
	double _abs;  /**< absolute value of the directed quantity */
};


/** @brief Dipole center
 *
 * Electrical dipole with dipole moment \f$\vec{p}\f$.
 */
class Dipole : public OrientedSite {
public:
	Dipole() : OrientedSite(0., 0., 0., 0., 0., 0., 0., 0.) {}
    /** @brief Constructor
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] eMyx     x coordinate of the dipole moments normal
     * \param[in] eMyy     y coordinate of the dipole moments normal
     * \param[in] eMyz     z coordinate of the dipole moments normal
     * \param[in] absMy     dipole moments absolute value
     */
    Dipole(double x, double y, double z, double eMyx, double eMyy, double eMyz, double absMy)
		: OrientedSite(x, y, z, 0., eMyx, eMyy, eMyz, absMy) {}

	/** @brief Read in XML configuration for a Dipole and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <site>
	     <!-- all Site class parameters -->
	     <dipolemoment>
	       <!-- either direct vector -->
	       <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z>
	       <!-- or polar angle represenation (DEG)-->
	       <theta>DOUBLE</theta> <phi>DOUBLE</phi>
	       <!-- -->
	       <abs>DOUBLE</abs> <!-- will overwrites abs obtained from an x y z specification -->
	     </dipolemoment>
	   </site>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) {
		Site::readXML(xmlconfig);
		Log::global_log->info() << "Site type: Dipole" << std::endl;
		bool bAngleInput = true;
		double theta, phi;
		bAngleInput = bAngleInput && xmlconfig.getNodeValueReduced("dipolemoment/theta", theta);
		bAngleInput = bAngleInput && xmlconfig.getNodeValueReduced("dipolemoment/phi",   phi);
		if(true == bAngleInput) {
			this->setOrientationVectorByPolarAngles(theta, phi);
		} else {
			xmlconfig.getNodeValueReduced("dipolemoment/x", _e[0]);
			xmlconfig.getNodeValueReduced("dipolemoment/y", _e[1]);
			xmlconfig.getNodeValueReduced("dipolemoment/z", _e[2]);
			normalize_e();
		}
		Log::global_log->info() << "Site parameters: dipole moment (ex, ey, ez) = (" << _e[0] << ", " << _e[1] << ", " << _e[2] << ")" << std::endl;
		xmlconfig.getNodeValueReduced("dipolemoment/abs", _abs);
		Log::global_log->info() << "Site parameters: dipole moment abs = " << _abs << std::endl;
	}

	double absMy() const { return abs(); }  /**< get the absolute value of the dipole moment. */

	/** set the absolute value of the dipole moment
	 * @note after setting components of the orientation vector normalization has to be performed
	 *       manually using the normalize_e() method.
	 */
	void setAbyMy(double my) { setAbs(my); }
};

/** @brief Quadrupole center
 *
 * Electrical quadrupole with quadrupole moment \f$\vec{Q}\f$.
 */
class Quadrupole : public OrientedSite {
public:
	Quadrupole() : OrientedSite(0., 0., 0., 0., 0., 0., 0., 0.) {}
    /** @brief Constructor
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] eQx      x coordinate of the quadrupole moments normal
     * \param[in] eQy      y coordinate of the quadrupole moments normal
     * \param[in] eQz      z coordinate of the quadrupole moments normal
     * \param[in] absQ     quadrupole moments absolute value
     */
	Quadrupole(double x, double y, double z, double eQx, double eQy, double eQz, double absQ)
		: OrientedSite(x, y, z, 0., eQx, eQy, eQz, absQ) {}

	/** @brief Read in XML configuration for a Dipole and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <site>
	     <!-- all Site class parameters -->
	     <quadrupolemoment>
	       <!-- either direct vector -->
	       <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z>
	       <!-- or polar angle represenation (DEG) -->
	       <theta>DOUBLE</theta> <phi>DOUBLE</phi>
	       <!-- -->
	       <abs>DOUBLE</abs> <!-- will overwrites abs obtained from an x y z specification -->
	     </quadrupolemoment>
	   </site>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) {
		Site::readXML(xmlconfig);
		Log::global_log->info() << "Site type: Quadrupole" << std::endl;
		bool bAngleInput = true;
		double theta, phi;
		bAngleInput = bAngleInput && xmlconfig.getNodeValueReduced("quadrupolemoment/theta", theta);
		bAngleInput = bAngleInput && xmlconfig.getNodeValueReduced("quadrupolemoment/phi",   phi);
		if(true == bAngleInput) {
			setOrientationVectorByPolarAngles(theta, phi);
		} else {
			xmlconfig.getNodeValueReduced("quadrupolemoment/x", _e[0]);
			xmlconfig.getNodeValueReduced("quadrupolemoment/y", _e[1]);
			xmlconfig.getNodeValueReduced("quadrupolemoment/z", _e[2]);
			normalize_e();
		}
		Log::global_log->info() << "Site parameters: quadrupole moment (ex, ey, ez) = (" << _e[0] << ", " << _e[1] << ", " << _e[2] << ")" << std::endl;
		xmlconfig.getNodeValueReduced("quadrupolemoment/abs", _abs);
		Log::global_log->info() << "Site parameters: quadrupole moment (abs) = " << _abs << std::endl;
	}

	double absQ() const { return abs(); }  /**< get the absolute value of the quadrupole moment. */

	/** set the absolute value of the quadrupole moment
	 * @note after setting components of the orientation vector normalization has to be performed
	 *       manually using the normalize_e() method.
	 */
	void setAbsQ(double q) { setAbs(q); }
};

class LJATMcenter : public LJcenter {
public:
	/** @brief Constructor */
	LJATMcenter(): _nu(0.) {}
	/** @brief Constructor
	 * \param[in] x		relative x coordinate
	 * \param[in] y		relative y coordinate
	 * \param[in] z		relative z coordinate
	 * \param[in] m		mass
	 * \param[in] nu	interaction strength
	 */
	LJATMcenter(double x, double y, double z,  double m, double epsilon, double sigma, double shift, double nu)
		: LJcenter(x, y, z, m, epsilon, sigma, shift), _nu(nu) {}

	/** @brief Read in XML configuration for an LJATMcenter and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <site type="LJATM">
	     <!-- all LJ126 Site class parameters -->
	     <nu>DOUBLE</nu>
	   </site>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) {
		LJcenter::readXML(xmlconfig);
		Log::global_log->info() << "Site type: Lennard-Jones combined with Axilrod-Teller-Muto" << std::endl;
		xmlconfig.getNodeValueReduced("nu", _nu);
		Log::global_log->info() << "Site parameters: nu = " << _nu << std::endl;
	}

	/// write to stream
	void write(std::ostream& ostrm) const {
		LJcenter::write(ostrm);
		ostrm << "\t" << nu();
	}

	/** Get the ATM interaction strength. */
	double nu() const { return _nu; }

	/** Set the ATM interaction strength. */
	void setNu(double nu) { _nu = nu; }

private:
	double _nu;
};
#endif  /* SITE_H_ */
