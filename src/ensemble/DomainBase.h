#ifndef DOMAIN_BASE_H
#define DOMAIN_BASE_H


class XMLfileUnits;

class DomainBase {
public:
	DomainBase(){};
	~DomainBase(){};

	virtual void readXML(XMLfileUnits& xmlconfig) = 0;

	/** Return volume of the domain */
	virtual double V() = 0;

	/** Return the lower corner of the domain.
	 *  Spanning a virtual box around the domain with sides parallel to the main axis x, y and z,
	 *  this is the upper corner of the minimal fictive box comprising the domain.
	 */
	double* rmin(){ return _rmin; }
	/** Return the upper corner of the domain.
	 *  Spanning a virtual box around the domain with sides parallel to the main axis x, y and z,
	 *  this is the lower corner of the minimal fictive box comprising the domain.
	 */
	double* rmax(){ return _rmax; }

	double length(int d) {
		return _rmax[d] - _rmin[d];
	}

protected:
	double _rmin[3]; /**< Lower corner of the minimal box around the domain. */
	double _rmax[3]; /**< Upper corner of the minimal box around the domain. */
};

#endif /* DOMAIN_BASE_H */
