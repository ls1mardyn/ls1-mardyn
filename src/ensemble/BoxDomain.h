#ifndef BOX_DOMAIN_H
#define BOX_DOMAIN_H

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
};

#endif /* BOX_DOMAIN_H */
