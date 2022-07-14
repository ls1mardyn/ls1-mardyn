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

	void setLength(int d, double l);
};

/* boxmin and boxmax not yet supported
* xml layout will be
* <boxMin> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </boxMin>
* <boxMax> <x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z> </boxMax>
* The default behavior of the domain is to set the boxMin coordinates to (0,0,0). So boxMin and boxMax are optional fields.
* If all 3 fields are given, boxMax will be ignored and calculated from l and boxMin. Otherwise the missing field will be calculated from the other two.
*/

#endif /* BOX_DOMAIN_H_ */
