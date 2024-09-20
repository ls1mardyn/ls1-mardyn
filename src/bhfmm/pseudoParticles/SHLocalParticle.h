/*
 * SHLocalParticle.h
 *
 *  Created on: Nov 27, 2014
 *      Author: tchipevn
 */

#ifndef SHLOCALPARTICLE_H_
#define SHLOCALPARTICLE_H_

#include "bhfmm/pseudoParticles/LocalParticle.h"
#include "bhfmm/expansions/SolidHarmonicsExpansion.h"
#include "bhfmm/utils/RotationParameter.h"
#include <map>

namespace bhfmm {

class SHLocalParticle: public LocalParticle {
public:
	SHLocalParticle(int order, bool initializeExpansionToZero = true);
	~SHLocalParticle();

	/**
	 * P2L operator
	 * @param position
	 * @param charge
	 */
	void addSource(const Vector3<double>& position, double charge);

	/**
	 * M2L operator
	 * @param multipole
	 * @param periodicShift - temporary addition for backwards compatibility with existing Periodic boundary conditions,
	 * @todo: remove, when new container for FMM is introduced
	 */
	void addMultipoleParticle(const MultipoleParticle& multipole, Vector3<double> periodicShift);

	/**
	 * rotation accelerated M2L operator
	 * @param multipole
	 * @param periodicShift - temporary addition for backwards compatibility with existing Periodic boundary conditions
	 * @param M2L_Wigner - lookup map for Wigner rotation matrices
	 */
	void addMultipoleParticle_Wigner(const MultipoleParticle& multipole, Vector3<double> periodicShift,
		double* cellWid,
		std::map<Vector3<int>, RotationParams, Vector3<int>::compare>& M2L_Wigner);

	/**
	 * L2L operator
	 * @param small
	 */
	void actOnLocalParticle(LocalParticle& small) const;

	/**
	 * rotation accelerated L2L operator
	 * @param small - the smaller particle to be added to the larger one (this)
	 * @param W_pos - Wigner matrix for positive rotation
	 * @param W_neg - Wigner matrix for negative (backward) rotation
	 * @param CosSinPhi - lookup for Cos(m*phi) and Sin(m*phi) values
	 * @param negate - reverse phi rotation (1 -> positive rotation, -1 -> negative rotation)
	 * @param magnitude - length of translation vector
	 */
	void actOnLocalParticle_Wigner(LocalParticle& small,
		const WignerMatrix& W_pos, const WignerMatrix& W_neg, const double* CosSinPhi, const int negate, const double& magnitude) const;

	/**
	 * L2P operator
	 * @param position
	 * @param charge
	 * @param potential stores resulting potential
	 * @param force stores resulting force
	 */
	void actOnTarget(const Vector3<double>& position, double charge, double& potential, Vector3<double>& force) const;

	void clear();
	int getNumEntries() const;


	const SolidHarmonicsExpansion& getConstExpansion() const {
		return _expansionM;
	}

	SolidHarmonicsExpansion& getExpansion() {
		return _expansionM;
	}

	void writeValuesToMPIBuffer(std::vector<double> buf, int& position) const {
		_expansionM.writeValuesToMPIBuffer(buf, position);
	}

	void readValuesFromMPIBuffer(std::vector<double> buf, int& position) {
		_expansionM.readValuesFromMPIBuffer(buf,position);
	}
	void addValuesFromMPIBuffer(std::vector<double> buf, int& position) {
		_expansionM.addValuesFromMPIBuffer(buf,position);
	}

private:
	SolidHarmonicsExpansion _expansionM;
};

} /* namespace bhfmm */

#endif /* SHLOCALPARTICLE_H_ */
