#ifndef SRC_UTILS_GENERATOR_MAXWELLVELOCITYASSIGNER_H_
#define SRC_UTILS_GENERATOR_MAXWELLVELOCITYASSIGNER_H_


#include <random>

#include "VelocityAssignerBase.h"

/** The MaxwellVelocityAssigner can be used to assign maxwell boltzmann distributed velocity vectors matching a given temperature.
 */
class MaxwellVelocityAssigner : public VelocityAssignerBase {
public:
	MaxwellVelocityAssigner(double T = 0, long seed = 0) : VelocityAssignerBase(T), _mt(seed), _normalDistribution(0.0, 1.0) {}
	~MaxwellVelocityAssigner() {}

	void assignVelocity(Molecule *molecule) {
		double v_abs = sqrt(/*kB=1*/ (1+molecule->component()->getRotationalDegreesOfFreedom()/3.)*T() / molecule->component()->m());
		double v[3];
		v[0] = v_abs * _normalDistribution(_mt);
		v[1] = v_abs * _normalDistribution(_mt);
		v[2] = v_abs * _normalDistribution(_mt);
		for(int d = 0; d < 3; d++) {
			molecule->setv(d, v[d]);
		}
	}
private:
	std::mt19937 _mt; //!< Mersenne twister used as input for the normal distribution
	std::normal_distribution<double> _normalDistribution;
};

#endif  // SRC_UTILS_GENERATOR_MAXWELLVELOCITYASSIGNER_H_
