/*
 * RDFDummyDecomposition.h
 *
 *  Created on: Jul 5, 2012
 *      Author: tijana
 */
#include "parallel/DomainDecompDummy.h"
#ifndef RDFDUMMYDECOMPOSITION_H_
#define RDFDUMMYDECOMPOSITION_H_



class RDFDummyDecomposition: public DomainDecompDummy {
public:
	RDFDummyDecomposition();
	virtual ~RDFDummyDecomposition();
	void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain);
	void exchangeMolecules(ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain);
};

#endif /* RDFDUMMYDECOMPOSITION_H_ */
