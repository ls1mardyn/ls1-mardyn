/*
 * SHMultipoleParticle.h
 *
 *  Created on: Nov 27, 2014
 *      Author: tchipevn
 */

#ifndef SHMULTIPOLEPARTICLE_H_
#define SHMULTIPOLEPARTICLE_H_

#include "bhfmm/pseudoParticles/MultipoleParticle.h"
#include "bhfmm/expansions/SolidHarmonicsExpansion.h"
#include "bhfmm/utils/WignerMatrix.h"

namespace bhfmm {

class SHMultipoleParticle: public MultipoleParticle {
public:
	SHMultipoleParticle(int order, bool initializeExpansionToZero = true);
	~SHMultipoleParticle();

	/**
	 * P2M operator
	 * @param position
	 * @param charge
	 */
	void addSource(const Vector3<double>& position, double charge);

	/**
	 * M2M operator
	 * @param small the smaller particle to be added to the larger one (this)
	 */
	void addMultipoleParticle(const MultipoleParticle& small);

	/**
	 * rotation accelerated M2M operator
	 * @param small - the smaller particle to be added to the larger one (this)
	 * @param W_pos - Wigner matrix for positive rotation
	 * @param W_neg - Wigner matrix for negative (backward) rotation
	 * @param CosSinPhi - lookup for Cos(m*phi) and Sin(m*phi) values
	 * @param negate - reverse phi rotation (1 -> positive rotation, -1 -> negative rotation)
	 * @param magnitude - length of translation vector
	 */
	void addMultipoleParticle_Wigner(const MultipoleParticle& small,
			const WignerMatrix& W_pos, const WignerMatrix& W_neg, const double* CosSinPhi, const int negate, const double& magnitude);

	/**
	 * M2L operator
	 * @param local
	 */
	void actOnLocalParticle(LocalParticle& local) const;

	/**
	 * M2P operator
	 * @param position
	 * @param charge
	 * @param potential stores resulting potential
	 * @param force stores resulting force
	 */
	void actOnTarget(const Vector3<double>& position, double charge, double& pot, Vector3<double>& force) const;

	void clear();
	int getNumEntries() const;

	const SolidHarmonicsExpansion& getConstExpansion() const {
		return _expansionL;
	}

	SolidHarmonicsExpansion& getExpansion() {
		return _expansionL;
	}

	void writeValuesToMPIBuffer(std::vector<double> buf, int& position) const {
		_expansionL.writeValuesToMPIBuffer(buf, position);
	}

	void readValuesFromMPIBuffer(std::vector<double> buf, int& position) {
		_expansionL.readValuesFromMPIBuffer(buf,position);
	}

	void addValuesFromMPIBuffer(std::vector<double> buf, int& position) {
		_expansionL.addValuesFromMPIBuffer(buf,position);
	}

private:
	SolidHarmonicsExpansion _expansionL;
};

} /* namespace bhfmm */

#endif /* SHMULTIPOLEPARTICLE_H_ */
