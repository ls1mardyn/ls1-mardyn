
#ifndef __LONGRANGECORRECTION_H__
#define __LONGRANGECORRECTION_H__

#include <cmath>

//class Domain;
//class Planar;
//class Homogeneous;
class XMLfileUnits;
class LongRangeCorrection{

public:
	LongRangeCorrection() {}
	virtual ~LongRangeCorrection() {}
//	void initializeLongRange();
	virtual void readXML(XMLfileUnits& xmlconfig) = 0;
	virtual void calculateLongRange() = 0;
/*
private:
	unsigned _type;
	Planar* _planar;
	Homogeneous* _homogen;
*/
  
};

#endif /* __LONGRANGECORRECTION_H__ */
