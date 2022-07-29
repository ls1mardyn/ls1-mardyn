#ifndef BOX_DOMAIN_H_
#define BOX_DOMAIN_H_

#include "DomainBase.h"

class XMLfileUnits;

class BoxDomain : public DomainBase {
public:
	BoxDomain();
	virtual ~BoxDomain(){}


	/** @brief Read in XML configuration for BoxDomain and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <domain type="box">
	     <lx>DOUBLE</lx>
	     <ly>DOUBLE</ly>
	     <lz>DOUBLE</lz>
	   </domain>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

	/** Return volume of the domain */
	virtual double V();

	virtual double length(int d);
	void setLength(int d, double l);
};

#endif /* BOX_DOMAIN_H_ */
