/*
 * FastMultipoleMethod.h
 *
 *  Created on: Feb 7, 2015
 *      Author: tchipev
 */

#ifndef FASTMULTIPOLEMETHOD_H_
#define FASTMULTIPOLEMETHOD_H_

#include "bhfmm/containers/PseudoParticleContainer.h"
#include "particleContainer/ParticleContainer.h"

namespace bhfmm {

class FastMultipoleMethod {
public:
	FastMultipoleMethod() : _order(-1), _LJCellSubdivisionFactor(0), _wellSeparated(0), _adaptive(false) {}
	~FastMultipoleMethod();

	/** @brief Read in XML configuration for FastMultipoleMethod and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <electrostatic type="FastMultipoleMethod">
	     <orderOfExpansions>UNSIGNED INTEGER</orderOfExpansions>
	     <LJCellSubdivisionFactor>INTEGER</LJCellSubdivisionFactor>
	   </electrostatic>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	void setParameters(unsigned LJSubdivisionFactor, int orderOfExpansions,
			bool periodic = true, bool adaptive = false);

	void init(double globalDomainLength[3], double bBoxMin[3],
			double bBoxMax[3], double LJCellLength[3], std::vector<int> neighbours);

	void computeElectrostatics(ParticleContainer * ljContainer);

	void printTimers();

private:
	int _order;
	unsigned _LJCellSubdivisionFactor;
	int _wellSeparated;
	int _adaptive;
	int _periodic;

	PseudoParticleContainer * _pseudoParticleContainer;

	VectorizedChargeP2PCellProcessor *_P2PProcessor;
	P2MCellProcessor *_P2MProcessor;
	L2PCellProcessor *_L2PProcessor;

};

} /* namespace bhfmm */

#endif /* FASTMULTIPOLEMETHOD_H_ */
