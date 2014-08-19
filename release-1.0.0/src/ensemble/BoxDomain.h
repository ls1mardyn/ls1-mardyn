#ifndef BOX_DOMAIN_H_
#define BOX_DOMAIN_H_

#include "DomainBase.h"

class XMLfileUnits;

class BoxDomain : public DomainBase {
public:
	BoxDomain();
	~BoxDomain(){}

	virtual void readXML(XMLfileUnits& xmlconfig);

	/** Return volume of the domain */
	virtual double V();

	virtual double length(int d);
	void setLength(int d, double l);
};

#endif /* BOX_DOMAIN_H_ */
