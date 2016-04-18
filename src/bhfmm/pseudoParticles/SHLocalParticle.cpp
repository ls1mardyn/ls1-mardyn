/*
 * SHLocalParticle.cpp
 *
 *  Created on: Nov 27, 2014
 *      Author: tchipevn
 */

#include "bhfmm/pseudoParticles/SHLocalParticle.h"
#include "bhfmm/pseudoParticles/SHMultipoleParticle.h"
#include <cassert>

namespace bhfmm {

SHLocalParticle::SHLocalParticle(int order, bool initializeExpansionToZero) :
		LocalParticle(), _expansionM(order, initializeExpansionToZero) {
	setOrder(order);
}

SHLocalParticle::~SHLocalParticle() {
}

void SHLocalParticle::addSource(const Vector3<double>& position, double charge) {
	// distance-vector FROM local particle TO point-source particle
	Vector3<double> r_pseudo_to_part = position - _center;

	// assert that the source does not lie in the local particle
	assert(r_pseudo_to_part.L2NormSquare() > _radiusSquared);

	_expansionM += charge * evaluateMOfR(_order, r_pseudo_to_part);
}

void SHLocalParticle::addMultipoleParticle(const MultipoleParticle& multipole, Vector3<double> periodicShift) {
	assert(multipole.getOrder() == _order);

	const SHMultipoleParticle& sh_multipole = static_cast<const SHMultipoleParticle&>(multipole);

	// compute periodically-shifted center
	Vector3<double> shifted_center = multipole.getCenter() + periodicShift;

	// distance-vector FROM local TO multipole
	Vector3<double> r_target_to_source = shifted_center - _center;
	//std::cout << r_target_to_source <<" \n";
	// assert that cells are not touching
	assert(r_target_to_source.L2Norm() >= _radius + multipole.getRadius());

	_expansionM += convoluteLM(setAtMinusR(sh_multipole.getConstExpansion()), evaluateMOfR(_order, r_target_to_source));
}

void SHLocalParticle::addMultipoleParticle_Wigner(const MultipoleParticle& multipole, Vector3<double> periodicShift,
		double* cellWid,
		std::map<Vector3<int>, RotationParams, Vector3<int>::compare>& M2L_Wigner) {
	assert(multipole.getOrder() == _order);

	const SHMultipoleParticle& sh_multipole = static_cast<const SHMultipoleParticle&>(multipole);

	// compute periodically-shifted center
	Vector3<double> shifted_center = multipole.getCenter() + periodicShift;

	// distance-vector FROM local TO multipole
	Vector3<double> r_target_to_source = shifted_center - _center;

	// assert that cells are not touching
	assert(r_target_to_source.L2Norm() >= _radius + multipole.getRadius());

	// distance-vector FROM big TO small along the z-axis
	Vector3<double> r_big_to_small_Z(0.0, 0.0, r_target_to_source.L2Norm());

	// compute index vector and add param to map //
	Vector3<int> idxVec(rint(r_target_to_source[0]/(cellWid[0])),
			rint(r_target_to_source[1]/(cellWid[1])),
			rint(r_target_to_source[2]/(cellWid[2])));

	// assert that Wigner matrices are present
	assert(M2L_Wigner.find(idxVec) != M2L_Wigner.end());

	RotationParams& param = M2L_Wigner[idxVec];

	SolidHarmonicsExpansion LM_trafo = rotatePhi(setAtMinusR(sh_multipole.getConstExpansion()), param.SinCos, 1); // pos. rotation
	LM_trafo = rotateThetaL(LM_trafo, param.W[0]);
	LM_trafo = convoluteLM_Z(LM_trafo, evaluateMOfR(_order, r_big_to_small_Z));
	LM_trafo = rotateThetaM(LM_trafo, param.W[1]);

	_expansionM += rotatePhi(LM_trafo, param.SinCos, -1); // neg. rotation
}

void SHLocalParticle::actOnLocalParticle(LocalParticle& small) const {
	assert(small.getOrder() == _order);

	SHLocalParticle& sh_small = static_cast<SHLocalParticle&>(small);

	// distance-vector FROM big local particle TO small local particle
	Vector3<double> r_big_to_small = small.getCenter() - _center;

	// assert that small is covered by big
	//assert(r_big_to_small.L2Norm() + small.getRadius() <= _radius);

	sh_small.getExpansion() += convoluteLM(evaluateLOfR(_order, r_big_to_small), _expansionM);
}

void SHLocalParticle::actOnLocalParticle_Wigner(LocalParticle& small,
		const WignerMatrix& W_pos, const WignerMatrix& W_neg, const double* CosSinPhi, const int negate, const double& magnitude) const {
	assert(small.getOrder() == _order);

	SHLocalParticle& sh_small = static_cast<SHLocalParticle&>(small);

	// distance-vector FROM big TO small along the z-axis
	Vector3<double> r_big_to_small_Z(0.0, 0.0, magnitude);

	// assert that small is covered by big
	assert(r_big_to_small_Z.L2Norm() + small.getRadius() <= _radius);

	SolidHarmonicsExpansion M_trafo = rotatePhi(_expansionM, CosSinPhi, negate); // pos. rotation
	M_trafo = rotateThetaM(M_trafo, W_pos);
	M_trafo = convoluteL_ZM(evaluateLOfR(_order, r_big_to_small_Z), M_trafo);
	M_trafo = rotateThetaM(M_trafo, W_neg);

	sh_small.getExpansion() += rotatePhi(M_trafo, CosSinPhi, -negate); // neg. rotation
}

void SHLocalParticle::actOnTarget(const Vector3<double>& position, double charge, double& pot,
		Vector3<double>& force) const {
	Vector3<double> r_pseudo_to_target = position - _center;

	// assert target does not lie in pseudo particle
	//assert(r_pseudo_to_target.L2NormSquare() >= _radiusSquared);

	SolidHarmonicsExpansion temp_L = evaluateLOfR(_order, r_pseudo_to_target);

	pot = potentialML(_expansionM, temp_L);
	force = charge * forceGradLAndM(temp_L, _expansionM);
}

void SHLocalParticle::clear() {
	_expansionM.clear();
}

int SHLocalParticle::getNumEntries() const {
	return _expansionM.getNumEntries();
}

} /* namespace bhfmm */

