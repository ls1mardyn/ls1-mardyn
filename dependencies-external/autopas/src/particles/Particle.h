/*
 * Particle.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLE_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLE_H_

namespace autopas {

class Particle {
public:
	Particle(double x = 0.0, double f = 0.0) :_x(x), _f(f) {}
	double _x;
	double _f;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLE_H_ */
