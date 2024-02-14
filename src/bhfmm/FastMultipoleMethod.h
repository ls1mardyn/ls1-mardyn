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
class UniformPseudoParticleContainer;
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
		currentMultipole, // or target
		sourceMultipole,
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
		P2PPreprocessSingleCell,
		P2PPostprocessSingleCell,
		P2Pc08StepBlock,
		P2MCompleteCell,
		M2MCompleteCell,
		M2LInitializeCell,
		M2LInitializeSource,
		M2LFinalizeCell,
		M2LTranslation,
		M2LPair2Way,
		L2LCompleteCell,
		L2PCompleteCell,
		Dummy
	};

private:
	int _order;
	unsigned _LJCellSubdivisionFactor;
	int _wellSeparated;
	bool _adaptive;
	bool _periodic;

	PseudoParticleContainer * _pseudoParticleContainer;

	VectorizedChargeP2PCellProcessor *_P2PProcessor;
	P2MCellProcessor *_P2MProcessor;
	L2PCellProcessor *_L2PProcessor;


#ifdef QUICKSCHED
	static void runner(int type, void *data);
	struct qsched *_scheduler;
#endif // QUICKSCEHD
};

} /* namespace bhfmm */

#endif /* FASTMULTIPOLEMETHOD_H_ */
