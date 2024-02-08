#ifndef COMPONENT_H_
#define COMPONENT_H_

#include <vector>
#include <string>

#include "molecules/Site.h"

/**
 * @brief Class implementing molecules as rigid rotators consisting out of different interaction sites (LJcenter, Charge, Dipole, Quadrupole).
 *
 * @author Martin Bernreuther
 */
class Component {
public:
	Component(unsigned int id = 0);

	/** @brief Read in XML configuration for a component and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <component id=INT name="STRING">
	     <site type="LJ126|Charge|Dipole|Quadrupole"> <!-- see Site class documentation --> </site>
	     <momentsofinertia> <Ixx>DOUBLE</Ixx> <Iyy>DOUBLE</Iyy> <Izz>DOUBLE</Izz> </momentsofinertia>
	   </component>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	/** set the component id */
	void setID(unsigned int id) { _id = id; }
	/** get component id */
	unsigned int ID() const { return _id; }

	/** get number of interaction sites */
	unsigned int numSites() const {
		return this->numLJcenters() + this->numCharges()
		                            + this->numDipoles()
		                            + this->numQuadrupoles();
	}

	/** get number of Lennard Jones interaction sites */
	unsigned int numLJcenters() const { return _ljcenters.size(); }
	/** get number of charge interaction sites */
	unsigned int numCharges() const { return _charges.size(); }
	/** get number of dipole interaction sites */
	unsigned int numDipoles() const { return _dipoles.size(); }
	/** get number of quadrupole interaction sites */
	unsigned int numQuadrupoles() const { return _quadrupoles.size(); }

	/** get mass of the molecule */
	double m() const { return _m; }
	double I11() const { return _Ipa[0]; }
	double I22() const { return _Ipa[1]; }
	double I33() const { return _Ipa[2]; }
	void setI11(double I) { _Ipa[0]=I; }
	void setI22(double I) { _Ipa[1]=I; }
	void setI33(double I) { _Ipa[2]=I; }

	/** get number of rotational degrees of freedom */
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

	void setNumMolecules(unsigned long num) { _numMolecules = num; }  /**< set the number of molecules for this component */
	void incNumMolecules() { ++_numMolecules; }  /**< increase the number of molecules for this component by 1 */
	void incNumMolecules(int N) { _numMolecules += N; }  /**< increase the number of molecules for this component by N */
	unsigned long getNumMolecules() const { return _numMolecules; }  /**< get the number of molecules (global) of the component */

	void addLJcenter(LJcenter& ljsite);
	void addLJcenter(
			double x, double y, double z, double m, double eps,
			double sigma, double rc = 0, bool TRUNCATED_SHIFTED = 0
	);
	void addCharge(Charge& chargesite);
	void addCharge(double x, double y, double z, double m, double q);
	void addDipole(Dipole& dipolesite);
	void addDipole(double x, double y, double z,
	               double eMyx, double eMyy, double eMyz, double eMyabs);
	void addQuadrupole(Quadrupole& quadrupolesite);
	void addQuadrupole(double x, double y, double z,
	                   double eQx, double eQy, double eQz, double eQabs);

	/** functions to fix the LJTS after rc is read in the simulation*/
	void updateAllLJcentersShift(double rc);

	/** delete the last site stored in the vector -- these are used by the external generators*/
	void deleteLJCenter() { _ljcenters.pop_back() ;}
	void deleteCharge() { _charges.pop_back() ;}
	void deleteDipole() { _dipoles.pop_back() ;}
	void deleteQuadrupole() { _quadrupoles.pop_back() ;}


	/**
	 * To be called after sites have been deleted or the properties of sites have been changed.
	 */
	void updateMassInertia();

	/** write information to stream */
	void write(std::ostream& ostrm) const;

	void writeVIM(std::ostream& ostrm);

	void setE_trans(double E) { _E_trans = E; }
	void setE_rot(double E) { _E_rot = E; }
	void setT(double T) { _T = T; }
	double E_trans() const { return _E_trans; }
	double E_rot() const { return _E_rot; }
	double E() const { return _E_trans + _E_rot; }
	double T() const { return _T; }
	void setName(std::string name) { _name = name; }
	std::string getName() const { return _name; }

	//! by Stefan Becker <stefan.becker@mv.uni-kl.de>
	//! needed by the MegaMol output format
	double getEps(unsigned int i) const {return _ljcenters[i].eps();}
	double getSigma(unsigned int i) const {return _ljcenters[i].sigma();}

	unsigned getLookUpId() const {
		return _lookUpID;
	}

	void setLookUpId(unsigned lookUpId) {
		_lookUpID = lookUpId;
	}

private:

	void updateMassInertia(Site& site);
	double calculateLJshift(double eps, double sigma, double rc) const;

	unsigned int _id; /**< component ID */
	// LJcenter,Dipole,Quadrupole have different size -> not suitable to store in a _Site_-array
	//std::vector<Site> _sites;
	// use separate vectors instead...
	std::vector<LJcenter> _ljcenters;
	std::vector<Charge> _charges;
	std::vector<Dipole> _dipoles;
	std::vector<Quadrupole> _quadrupoles;
	/* for performance reasons better(?) omit Site-class indirection and use cached values */
	double _m; /**< total mass */
	/** moments of inertia tensor
	 * \f$ I_{xx}, I_{yy}, I_{zz}, I_{xy}, I_{xz}, I_{yz} \f$
	 */
	double _I[6];
	double _Ipa[3]; /**< moments of inertia for principal axes */
	unsigned long _rot_dof; /**< number of rotational degrees of freedom */

	/* cached values, set by ensemble class! */
	unsigned long _numMolecules; // number of molecules for this molecule type
	double _E_trans; // translational energy
	double _E_rot; // rotational energy
	double _T; // temperature

	bool _isStockmayer; //Checks whether component is a Stockmayer fluid to determine moments of inertia

	std::string _name; /**< name of the component/molecule type */

	/**
	 * for use by the Vectorization:
	 * a look-up table is set up there,
	 * which needs to know about the sites of all present components
	 *
	 * \brief One LJ center enumeration start index for each component.
	 * \details All the LJ centers of all components are enumerated.<br>
	 * Comp1 gets indices 0 through n1 - 1, Comp2 n1 through n2 - 1 and so on.<br>
	 * This is necessary for finding the respective parameters for each interaction<br>
	 * between two centers.
	 */
	unsigned _lookUpID;
};

std::ostream& operator<<(std::ostream& stream, const Component& component);

#endif /* COMPONENT_H_ */
