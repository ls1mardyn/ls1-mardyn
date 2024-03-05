
#ifndef __LONGRANGECORRECTION_H__
#define __LONGRANGECORRECTION_H__

#include <cmath>

class Domain;
class DomainDecompBase;
//class Planar;
//class Homogeneous;
class XMLfileUnits;
class LongRangeCorrection{

public:
	LongRangeCorrection() {}
	virtual ~LongRangeCorrection() {}
	virtual void init() = 0;
	virtual void readXML(XMLfileUnits& xmlconfig) = 0;
	virtual void calculateLongRange() = 0;
	virtual void writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) = 0;
/*
private:
	unsigned _type;
	Planar* _planar;
	Homogeneous* _homogen;
*/

};

#endif /* __LONGRANGECORRECTION_H__ */
