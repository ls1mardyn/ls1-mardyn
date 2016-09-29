/*
 * PseudoParticle.h
 *
 *  Created on: Nov 27, 2014
 *      Author: tchipevn
 */

#ifndef PSEUDOPARTICLE_H_
#define PSEUDOPARTICLE_H_

#include "bhfmm/utils/Vector3.h"

namespace bhfmm {

/**
 * Base class for multipole and local pseudo-particles. Define
 * center, radius and order. Radius will be more extensively used (and tested) in the adaptive case.
 */
class PseudoParticle {
public:
	PseudoParticle() :
			_center(0.0), _radius(0.0), _radiusSquared(0.0), _order(0) {
	}

	virtual ~PseudoParticle() {
	}

	/**
	 * @return the center of the particle
	 */
	const Vector3<double>& getCenter() const {
		return _center;
	}

	/**
	 * set the center of the particle
	 * @param center
	 */
	void setCenter(const Vector3<double>& center) {
		_center = center;
	}

	/**
	 * @todo: is this needed?
	 * @return the order of the particle
	 */
	int getOrder() const {
		return _order;
	}

	/**
	 * set the order
	 * @param order
	 */
	void setOrder(int order) {
		_order = order;
	}

	/**
	 * get the radius of the pseudo particle
	 * @return radius
	 */
	double getRadius() const {
		return _radius;
	}

	/**
	 * set the radius
	 * @param radius
	 */
	void setRadius(double radius) {
		_radius = radius;
		_radiusSquared = radius * radius;
	}

	/**
	 * set expansions to zero for a new iteration
	 */
	virtual void clear() = 0;

	virtual int getNumEntries() const = 0;

protected:
	Vector3<double> _center;

	double _radius;

	double _radiusSquared;

	// maybe not needed:
	int _order;

};

} /* namespace bhfmm */

#endif /* PSEUDOPARTICLE_H_ */
