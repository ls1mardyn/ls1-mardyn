#ifndef SRC_UTILS_GENERATOR_VELOCITYASSIGNERBASE_H_
#define SRC_UTILS_GENERATOR_VELOCITYASSIGNERBASE_H_

#include "molecules/Molecule.h"

/** The VelocityAssignerBase implements the gernal functionality and interface to assign velocity vectors mathing to a given temperature.
 */
class VelocityAssignerBase {
public:
	VelocityAssignerBase(double T = 0) : _T(T) {}
	virtual ~VelocityAssignerBase(){}
	void setTemperature(double T) { _T = T; }
	double T() { return _T; }
	virtual void assignVelocity(Molecule *molecule) = 0;
private:
	double _T;  //!< coressponding target temperature
};

#endif  // SRC_UTILS_GENERATOR_VELOCITYASSIGNERBASE_H_
