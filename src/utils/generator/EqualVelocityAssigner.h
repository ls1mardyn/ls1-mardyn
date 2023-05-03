#ifndef SRC_UTILS_GENERATOR_EQUALVELOCITYASSIGNER_H_
#define SRC_UTILS_GENERATOR_EQUALVELOCITYASSIGNER_H_

#include <random>

#include "VelocityAssignerBase.h"

/** The VelocityAssigner can be used to assign normally distributed velocity vectors with absolute value matching a given temperature.
 */
class EqualVelocityAssigner : public VelocityAssignerBase {
public:
	EqualVelocityAssigner(double T = 0, long seed = 0) : VelocityAssignerBase(T), _mt(seed), _uniformDistribution(0, 1) {}
	~EqualVelocityAssigner(){}

	void assignVelocity(Molecule *molecule) {
		double v_abs = sqrt(/*kB=1*/ (3+molecule->component()->getRotationalDegreesOfFreedom())*T() / molecule->component()->m());
		/* pick angels for uniform distributino on S^2. */
		double phi, theta;
		phi   = 2*M_PI * _uniformDistribution(_mt);
		theta = acos(2 * _uniformDistribution(_mt) - 1);
		double v[3];
		v[0] = v_abs * sin(phi);
		v[1] = v_abs * cos(phi) * sin(theta);
		v[2] = v_abs * cos(phi) * cos(theta);
		for(int d = 0; d < 3; d++) {
			molecule->setv(d, v[d]);
		}
	}
private:
	std::mt19937 _mt; //!< Mersenne twister used as input for the uniform distribution
	std::uniform_real_distribution<double> _uniformDistribution;
};

#endif  // SRC_UTILS_GENERATOR_EQUALVELOCITYASSIGNER_H_
