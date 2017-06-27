/*
 * FastMultipoleMethod.cpp
 *
 *  Created on: Feb 7, 2015
 *      Author: tchipev
 */

#include <utils/threeDimensionalMapping.h>
#include "FastMultipoleMethod.h"
#include "Simulation.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "bhfmm/containers/UniformPseudoParticleContainer.h"
#include "bhfmm/containers/AdaptivePseudoParticleContainer.h"
#include "utils/xmlfileUnits.h"
#ifdef USE_VT
#include "VT.h"
#endif

using Log::global_log;
using std::endl;

namespace bhfmm {

FastMultipoleMethod::~FastMultipoleMethod() {
	delete _pseudoParticleContainer;
	delete _P2PProcessor;
	delete _P2MProcessor;
	delete _L2PProcessor;
#ifdef QUICKSCHED
    qsched_free(_scheduler);
    delete(_scheduler);
#endif
}

void FastMultipoleMethod::readXML(XMLfileUnits& xmlconfig) {

	xmlconfig.getNodeValue("orderOfExpansions", _order);
	global_log->info() << "FastMultipoleMethod: orderOfExpansions: " << _order << endl;

	xmlconfig.getNodeValue("LJCellSubdivisionFactor", _LJCellSubdivisionFactor);
	global_log->info() << "FastMultipoleMethod: LJCellSubdivisionFactor: " << _LJCellSubdivisionFactor << endl;

	xmlconfig.getNodeValue("adaptiveContainer", _adaptive);
	if (_adaptive == 1) {
		global_log->warning() << "FastMultipoleMethod: adaptiveContainer is not debugged yet and certainly delivers WRONG results!" << endl;
		global_log->warning() << "Unless you are in the process of debugging this container, please stop the simulation and restart with the uniform one" << endl;
	} else {
		global_log->info() << "FastMultipoleMethod: UniformPseudoParticleSelected " << endl;
	}

	xmlconfig.getNodeValue("systemIsPeriodic", _periodic);
	if (_periodic == 0) {
		global_log->warning() << "FastMultipoleMethod: periodicity is turned off!" << endl;
	} else {
		global_log->info() << "FastMultipoleMethod: Periodicity is on." << endl;
	}
}

void FastMultipoleMethod::setParameters(unsigned LJSubdivisionFactor,
		int orderOfExpansions, bool periodic, bool adaptive) {
	_LJCellSubdivisionFactor = LJSubdivisionFactor;
	_order = orderOfExpansions;
	_periodic = periodic;
	_adaptive = adaptive;
}

void FastMultipoleMethod::init(double globalDomainLength[3], double bBoxMin[3],
		double bBoxMax[3], double LJCellLength[3], ParticleContainer* ljContainer) {

	if (_LJCellSubdivisionFactor != 1
			and _LJCellSubdivisionFactor != 2
			and _LJCellSubdivisionFactor != 4
			and _LJCellSubdivisionFactor != 8) {
		global_log->error() << "Fast Multipole Method: bad subdivision factor:"
				<< _LJCellSubdivisionFactor << endl;
		global_log->error() << "expected 1,2,4 or 8" << endl;
		Simulation::exit(5);
	}
	global_log->info()
			<< "Fast Multipole Method: each LJ cell will be subdivided in "
			<< pow(_LJCellSubdivisionFactor, 3)
			<< " cells for electrostatic calculations in FMM" << endl;

	_P2PProcessor = new VectorizedChargeP2PCellProcessor(
			*(global_simulation->getDomain()));
#ifdef QUICKSCHED
    _scheduler = new struct qsched;
    qsched_init(_scheduler, mardyn_get_max_threads(), qsched_flag_none);
#endif // QUICKSCEHD
	if (not _adaptive) {
		_pseudoParticleContainer = new UniformPseudoParticleContainer(globalDomainLength,
                                                                      bBoxMin,
                                                                      bBoxMax,
                                                                      LJCellLength,
                                                                      _LJCellSubdivisionFactor,
                                                                      _order,
                                                                      ljContainer,
                                                                      _periodic
#ifdef QUICKSCHED
                                                                    , _scheduler
#endif
        );

	} else {
		// TODO: Debugging in Progress!
#if defined(ENABLE_MPI)
		global_log->error() << "MPI in combination with adaptive is not supported yet" << endl;
		Simulation::exit(-1);
#endif
		//int threshold = 100;
		_pseudoParticleContainer = new AdaptivePseudoParticleContainer(
				globalDomainLength, _order, LJCellLength,
				_LJCellSubdivisionFactor, _periodic);
	}

	_P2MProcessor = new P2MCellProcessor(_pseudoParticleContainer);
	_L2PProcessor = new L2PCellProcessor(_pseudoParticleContainer);

    contextFMM = this;
}

void FastMultipoleMethod::computeElectrostatics(ParticleContainer* ljContainer) {
#ifdef USE_VT
	VT_traceon();
#endif
	// build
	_pseudoParticleContainer->build(ljContainer);

	// clear expansions
	_pseudoParticleContainer->clear();

	// P2M, M2P
	_pseudoParticleContainer->upwardPass(_P2MProcessor);

	// M2L, P2P
	if (_adaptive) {
		_P2PProcessor->initTraversal();
	}
#ifdef QUICKSCHED
    qsched_run(_scheduler, mardyn_get_max_threads(), runner);
#endif
	_pseudoParticleContainer->horizontalPass(_P2PProcessor);

	// L2L, L2P
	_pseudoParticleContainer->downwardPass(_L2PProcessor);
#ifdef USE_VT
	VT_traceoff();
#endif
}

void FastMultipoleMethod::printTimers() {
	_P2PProcessor->printTimers();
	_P2MProcessor->printTimers();
	_L2PProcessor->printTimers();
	_pseudoParticleContainer->printTimers();
	//global_simulation->printTimers("CELL_PROCESSORS");
	//global_simulation->printTimers("UNIFORM_PSEUDO_PARTICLE_CONTAINER");
}

#ifdef  QUICKSCHED
void FastMultipoleMethod::runner(int type, void *data) {
	switch (type) {
		case PreprocessCell:{
            ParticleCellPointers * cell = ((ParticleCellPointers **)data)[0];
            contextFMM->_P2PProcessor->preprocessCell(*cell);
			break;
		} /* PreprocessCell */
		case PostprocessCell:{
            ParticleCellPointers * cell = ((ParticleCellPointers **)data)[0];
            contextFMM->_P2PProcessor->postprocessCell(*cell);
			break;
		} /* PostprocessCell */
		case P2P:{
            // TODO optimize calculation order (1. corners 2. edges 3. rest) and gradually release resources
            unsigned long x = ((unsigned long *) data)[0];
            unsigned long y = ((unsigned long *) data)[1];
            unsigned long z = ((unsigned long *) data)[2];
            unsigned long blockX = ((unsigned long *) data)[3];
            unsigned long blockY = ((unsigned long *) data)[4];
            unsigned long blockZ = ((unsigned long *) data)[5];
            LeafNodesContainer *contextContainer = ((LeafNodesContainer **) data)[6];

            // traverse over block
            for (unsigned long i = 0; i < blockX - 1
                                      && i < contextContainer->getNumCellsPerDimension()[0] - 1; ++i) {
                for (unsigned long j = 0; j < blockY - 1
                                          && j < contextContainer->getNumCellsPerDimension()[1] - 1; ++j) {
                    for (unsigned long k = 0; k < blockZ - 1
                                              && k < contextContainer->getNumCellsPerDimension()[2] - 1; ++k) {

                        //process cell
                        long baseIndex = contextContainer->cellIndexOf3DIndex(x + i,
                                                                     y + j,
                                                                     z + k);
                        contextFMM->_P2PProcessor->processCell(contextContainer->getCells()[baseIndex]);
                    }
                }
            }
			break;
		} /* P2P */
		case Dummy:{
            // do nothing, only serves for synchronization
			break;
		} /* P2P */
		default:
			global_log->error() << "Undefined Quicksched task type: " << type << std::endl;
	}
}
#endif /* QUICKSCEHD */
} /* namespace bhfmm */

