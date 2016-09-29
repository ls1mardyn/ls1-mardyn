/*
 * SHMultipoleParticle.cpp
 *
 *  Created on: Nov 27, 2014
 *      Author: tchipevn
 */

#include "bhfmm/pseudoParticles/SHMultipoleParticle.h"
#include "bhfmm/pseudoParticles/SHLocalParticle.h"
#include <cassert>

namespace bhfmm {

SHMultipoleParticle::SHMultipoleParticle(int order, bool initializeExpansionToZero) :
		MultipoleParticle(), _expansionL(order, initializeExpansionToZero) {
	setOrder(order);
}

SHMultipoleParticle::~SHMultipoleParticle() {
}

void SHMultipoleParticle::addSource(const Vector3<double>& position, double charge) {
	// distance-vector FROM multipole particle TO point-source particle
	Vector3<double> r_pseudo_to_part = position - _center;

	/* ToDo: fails: corner case centers should go to different cells */
	//assert(r_pseudo_to_part.L2NormSquare() < _radiusSquared);

	_expansionL += charge * evaluateLOfR(_order, r_pseudo_to_part);
}

void SHMultipoleParticle::addMultipoleParticle(const MultipoleParticle& small) {
	assert(small.getOrder() == _order);

	const SHMultipoleParticle& sh_small = static_cast<const SHMultipoleParticle&>(small);

	// distance-vector FROM big TO small
	Vector3<double> r_big_to_small = small.getCenter() - _center;

	// assert that small is covered by big
	/* ToDo: fails: maybe this can be covered by adding a small eps: _radius + eps */
	//assert(r_big_to_small.L2Norm() + small.getRadius() <= _radius);

	_expansionL += convoluteLL(sh_small._expansionL, evaluateLOfR(_order, r_big_to_small));

}

void SHMultipoleParticle::addMultipoleParticle_Wigner(const MultipoleParticle& small,
		const WignerMatrix& W_pos, const WignerMatrix& W_neg, const double* CosSinPhi, const int negate, const double& magnitude) {
	assert(small.getOrder() == _order);

	const SHMultipoleParticle& sh_small = static_cast<const SHMultipoleParticle&>(small);

	// distance-vector FROM big TO small along the z-axis
	Vector3<double> r_big_to_small_Z(0.0, 0.0, magnitude);

	// assert that small is covered by big
	assert(r_big_to_small_Z.L2Norm() + small.getRadius() <= _radius);

	SolidHarmonicsExpansion L_trafo = rotatePhi(sh_small._expansionL, CosSinPhi, negate); // pos. rotation
	L_trafo = rotateThetaL(L_trafo, W_pos);
	L_trafo = convoluteLL_Z(L_trafo, evaluateLOfR(_order, r_big_to_small_Z));
	L_trafo = rotateThetaL(L_trafo, W_neg);

	_expansionL += rotatePhi(L_trafo, CosSinPhi, -negate); // neg. rotation
}

void SHMultipoleParticle::actOnLocalParticle(LocalParticle& local) const {
	assert(local.getOrder() == _order);

	SHLocalParticle& sh_local = static_cast<SHLocalParticle&>(local);

	// distance-vector FROM local TO source
	Vector3<double> r_target_to_source = _center - local.getCenter();

	// assert that cells are not touching
	assert(r_target_to_source.L2Norm() >= _radius + local.getRadius());

	sh_local.getExpansion() += convoluteLM(setAtMinusR(_expansionL), evaluateMOfR(_order, r_target_to_source));
}

void SHMultipoleParticle::actOnTarget(const Vector3<double>& position, double /*charge*/, double& pot,
		Vector3<double>& force) const {
	Vector3<double> r_pseudo_to_target = position - _center;

	// assert that point-target particle does not lie in pseudo particle
	assert(r_pseudo_to_target.L2NormSquare() >= _radiusSquared);

	// for the force, generate expansion of order (ord+1),
	// in order to fully utilise the current L expansion (due to particular stencil of gradient)
	SolidHarmonicsExpansion temp_M = evaluateMOfR(_order + 1, r_pseudo_to_target);

	pot = potentialML(temp_M, _expansionL);
	force = forceLAndGradM(_expansionL, temp_M);
}

void SHMultipoleParticle::clear() {
	_expansionL.clear();
}

int SHMultipoleParticle::getNumEntries() const {
	return _expansionL.getNumEntries();
}

} /* namespace bhfmm */

