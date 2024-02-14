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
	Log::global_log->info() << "FastMultipoleMethod: orderOfExpansions: " << _order << std::endl;

	xmlconfig.getNodeValue("LJCellSubdivisionFactor", _LJCellSubdivisionFactor);
	Log::global_log->info() << "FastMultipoleMethod: LJCellSubdivisionFactor: " << _LJCellSubdivisionFactor << std::endl;

	xmlconfig.getNodeValue("adaptiveContainer", _adaptive);
	if (_adaptive) {
		Log::global_log->warning() << "FastMultipoleMethod: adaptiveContainer is not debugged yet and certainly delivers WRONG results!" << std::endl;
		Log::global_log->warning() << "Unless you are in the process of debugging this container, please stop the simulation and restart with the uniform one" << std::endl;
	} else {
		Log::global_log->info() << "FastMultipoleMethod: UniformPseudoParticleSelected " << std::endl;
	}

	xmlconfig.getNodeValue("systemIsPeriodic", _periodic);
	if (_periodic) {
		Log::global_log->info() << "FastMultipoleMethod: Periodicity is on." << std::endl;
	} else {
		Log::global_log->warning() << "FastMultipoleMethod: periodicity is turned off!" << std::endl;
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
		Log::global_log->error() << "Fast Multipole Method: bad subdivision factor:"
				<< _LJCellSubdivisionFactor << std::endl;
		Log::global_log->error() << "expected 1,2,4 or 8" << std::endl;
		Simulation::exit(5);
	}
	Log::global_log->info()
			<< "Fast Multipole Method: each LJ cell will be subdivided in "
			<< pow(_LJCellSubdivisionFactor, 3)
			<< " cells for electrostatic calculations in FMM" << std::endl;

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
#ifdef TASKTIMINGPROFILE
#ifdef QUICKSCHED
        global_simulation->getTaskTimingProfiler()->init(_scheduler->count);
#else
        Log::global_log->warning() << "Profiling tasks without Quicksched not implemented!" << std::endl;
#endif
#endif

	} else {
		// TODO: Debugging in Progress!
#if defined(ENABLE_MPI)
		Log::global_log->error() << "MPI in combination with adaptive is not supported yet" << std::endl;
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
	// build
	_pseudoParticleContainer->build(ljContainer);

	// clear expansions
	_pseudoParticleContainer->clear();

#ifdef QUICKSCHED
	_P2MProcessor->initTraversal();
	_P2PProcessor->initTraversal();
	_L2PProcessor->initTraversal();

	global_simulation->timers()->start(("QUICKSCHED"));
	qsched_run(_scheduler, mardyn_get_max_threads(), runner);
	global_simulation->timers()->stop(("QUICKSCHED"));

	_L2PProcessor->endTraversal();
	_P2PProcessor->endTraversal();
	_P2MProcessor->endTraversal();

	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_FMM_COMPLETE");
#else
    // P2M, M2P
	_pseudoParticleContainer->upwardPass(_P2MProcessor);
	if (_adaptive) {
		_P2PProcessor->initTraversal();
	}
	// M2L, P2P
	_pseudoParticleContainer->horizontalPass(_P2PProcessor);
	// L2L, L2P
	_pseudoParticleContainer->downwardPass(_L2PProcessor);
#endif

}

void FastMultipoleMethod::printTimers() {
	_P2PProcessor->printTimers();
	_P2MProcessor->printTimers();
	_L2PProcessor->printTimers();
	_pseudoParticleContainer->printTimers();
	global_simulation->timers()->printTimers("CELL_PROCESSORS");
	global_simulation->timers()->printTimers("UNIFORM_PSEUDO_PARTICLE_CONTAINER");
}

#ifdef QUICKSCHED
void FastMultipoleMethod::runner(int type, void *data) {
#ifdef FMM_FFT
	struct qsched_payload *payload = (qsched_payload *)data;
#ifdef TASKTIMINGPROFILE
    auto startTime = global_simulation->getTaskTimingProfiler()->start();
#endif /* TASKTIMINGPROFILE */
	switch (type) {
		case P2PPreprocessSingleCell:{
            contextFMM->_P2PProcessor->preprocessCell(*payload->cell.pointer);
			break;
		} /* PreprocessCell */
		case P2PPostprocessSingleCell:{
            contextFMM->_P2PProcessor->postprocessCell(*payload->cell.pointer);
			break;
		} /* PostprocessCell */
		case P2Pc08StepBlock:{
            // TODO optimize calculation order (1. corners 2. edges 3. rest) and gradually release resources
			int                x                 = payload->cell.coordinates[0];
			int                y                 = payload->cell.coordinates[1];
			int                z                 = payload->cell.coordinates[2];
			LeafNodesContainer *contextContainer = payload->leafNodesContainer;

			long baseIndex = contextContainer->cellIndexOf3DIndex(x, y, z);
			contextContainer->c08Step(baseIndex, *contextFMM->_P2PProcessor);

			break;
		} /* P2Pc08StepBlock */
		// used for scheme Pair2Way
        case M2LInitializeCell: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            if (contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole].occ == 0)
                break;

            double radius = contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .local
                    .getRadius();

            FFTAccelerableExpansion &source = static_cast<bhfmm::SHMultipoleParticle &>(contextContainer
                    ->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .multipole)
                    .getExpansion();
            FFTAccelerableExpansion &target = static_cast<bhfmm::SHLocalParticle &>(contextContainer
                    ->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .local)
                    .getExpansion();
            contextContainer->getFFTAcceleration()->FFT_initialize_Source(source, radius);
            contextContainer->getFFTAcceleration()->FFT_initialize_Target(target);
			break;
        } /* M2LInitializeCell */
		// used for scheme CompleteTarget
        case M2LInitializeSource: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            if (contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole].occ == 0)
                break;

            double radius = contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .local
                    .getRadius();

            FFTAccelerableExpansion &source = static_cast<bhfmm::SHMultipoleParticle &>(contextContainer
                    ->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .multipole)
                    .getExpansion();
            contextContainer->getFFTAcceleration()->FFT_initialize_Source(source, radius);
            break;
        } /* M2LInitializeSource */
		// used for scheme Pair2Way
        case M2LFinalizeCell: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            if (contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole].occ == 0)
                break;

            double radius = contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .local
                    .getRadius();

            FFTAccelerableExpansion &target = static_cast<bhfmm::SHLocalParticle &>(contextContainer
                    ->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .local)
                    .getExpansion();
            contextContainer->getFFTAcceleration()->FFT_finalize_Target(target, radius);
            break;
        } /* M2LFinalizeCell */
		// used for scheme CompleteTarget
        case M2LTranslation: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            contextContainer->M2LCompleteCell(payload->currentMultipole,
                                              payload->currentLevel,
                                              payload->currentEdgeLength);
			break;
		} /* M2LTranslation */
		// used for scheme Pair2Way
		case M2LPair2Way: {
			UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

			contextContainer->M2LPair2Way(payload->sourceMultipole,
										  payload->currentMultipole,
										  payload->currentLevel,
										  payload->currentEdgeLength);
			break;
		}
        case L2LCompleteCell: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            contextContainer->L2LCompleteCell(payload->sourceMultipole,
											  payload->currentLevel,
											  payload->currentEdgeLength);
			break;
        }
		case L2PCompleteCell: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            contextContainer->L2PCompleteCell(payload->currentMultipole);
            break;
		}
        case M2MCompleteCell: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            contextContainer->M2MCompleteCell(payload->currentMultipole,
                                              payload->currentLevel,
                                              payload->currentEdgeLength);
            break;
        }
        case P2MCompleteCell: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            contextContainer->P2MCompleteCell(payload->currentMultipole);
            break;
        }
		case Dummy: {
			// do nothing, only serves for synchronization
			break;
		} /* Dummy */
		default:
			Log::global_log->error() << "Undefined Quicksched task type: " << type << std::endl;
    }
#ifdef TASKTIMINGPROFILE
    global_simulation->getTaskTimingProfiler()->stop(startTime, type);
#endif /* TASKTIMINGPROFILE */
#else
#pragma omp critical
	{
	Log::global_log->error() << "Quicksched runner without FMM_FFT not implemented!" << std::endl;
	Simulation::exit(1);
	}
#endif /* FMM_FFT */
}
#endif /* QUICKSCEHD */
} /* namespace bhfmm */

