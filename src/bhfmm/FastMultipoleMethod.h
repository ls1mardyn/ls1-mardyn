/*
 * FastMultipoleMethod.h
 *
 *  Created on: Feb 7, 2015
 *      Author: tchipev
 */

#ifndef FASTMULTIPOLEMETHOD_H_
#define FASTMULTIPOLEMETHOD_H_

#include <bhfmm/containers/LeafNodesContainer.h>
#include <bhfmm/containers/UniformPseudoParticleContainer.h>
#include "bhfmm/containers/PseudoParticleContainer.h"
#include "particleContainer/ParticleContainer.h"
#ifdef QUICKSCHED
#include "quicksched.h"
#endif

namespace bhfmm {
class FastMultipoleMethod;
// needed for static runner()
static FastMultipoleMethod *contextFMM;

#ifdef QUICKSCHED
    //TODO more unions to save space
struct qsched_payload {
	// Stuff for P2P and pre/post processing
    union {
        ParticleCellPointers *pointer;
	    int coordinates[3];
    } cell;
	int taskBlockSize[3];
	LeafNodesContainer *leafNodesContainer;

	// Stuff for M2L and init/finalize FFT

	int currentLevel,
		currentMultipole,
        currentEdgeLength;
	UniformPseudoParticleContainer *uniformPseudoParticleContainer;
};
#endif
class FastMultipoleMethod {
public:
	FastMultipoleMethod() : _order(-1),
                            _LJCellSubdivisionFactor(0),
                            _wellSeparated(0),
                            _adaptive(false)
    {}
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
			double bBoxMax[3], double LJCellLength[3], ParticleContainer* ljContainer);

	void computeElectrostatics(ParticleContainer * ljContainer);

	void printTimers();

	enum taskType {
		Dummy, P2PPreprocessSingleCell, P2PPostprocessSingleCell, P2Pc08StepBlock, M2LInitializePair, M2LFinalizePair, M2LTranslation
	};
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


#ifdef QUICKSCHED
	static void runner(int type, void *data);
	struct qsched *_scheduler_p2p;
	struct qsched *_scheduler_m2l;
#endif // QUICKSCEHD
};

} /* namespace bhfmm */

#endif /* FASTMULTIPOLEMETHOD_H_ */
