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
    qsched_free(_scheduler_p2p);
    delete(_scheduler_p2p);
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
    _scheduler_p2p = new struct qsched;
    _scheduler_m2l = new struct qsched;
    qsched_init(_scheduler_p2p, mardyn_get_max_threads(), qsched_flag_none);
    qsched_init(_scheduler_m2l, mardyn_get_max_threads(), qsched_flag_none);
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
                                                                    , _scheduler_p2p
                                                                    , _scheduler_m2l
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
    qsched_run(_scheduler_p2p, mardyn_get_max_threads(), runner);
    qsched_run(_scheduler_m2l, mardyn_get_max_threads(), runner);
#else
	_pseudoParticleContainer->horizontalPass(_P2PProcessor);
#endif

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
	struct qsched_payload *payload = (qsched_payload *)data;
	switch (type) {
		case PreprocessCell:{
            contextFMM->_P2PProcessor->preprocessCell(*payload->cell.pointer);
			break;
		} /* PreprocessCell */
		case PostprocessCell:{
            contextFMM->_P2PProcessor->postprocessCell(*payload->cell.pointer);
			break;
		} /* PostprocessCell */
		case P2P:{
            // TODO optimize calculation order (1. corners 2. edges 3. rest) and gradually release resources
			int                x                 = payload->cell.coordinates[0];
			int                y                 = payload->cell.coordinates[1];
			int                z                 = payload->cell.coordinates[2];
			LeafNodesContainer *contextContainer = payload->leafNodesContainer;

			long baseIndex = contextContainer->cellIndexOf3DIndex(x, y, z);
			contextContainer->c08Step(baseIndex, *contextFMM->_P2PProcessor);

			break;
		} /* P2P */
        case FFTInitialize: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            if (contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole].occ == 0)
                return;

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
        } /* FFTInitialize */
        case FFTFinalize: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            if (contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole].occ == 0)
                return;

            double radius = contextContainer->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .local
                    .getRadius();

            FFTAccelerableExpansion &target = static_cast<bhfmm::SHLocalParticle &>(contextContainer
                    ->getMpCellGlobalTop()[payload->currentLevel][payload->currentMultipole]
                    .local)
                    .getExpansion();
            contextContainer->getFFTAcceleration()->FFT_finalize_Target(target, radius);
            break;
        } /* FFTFinalize */
        case M2LFourier: {
            UniformPseudoParticleContainer *contextContainer = payload->uniformPseudoParticleContainer;

            contextContainer->M2LStep<true, true, false, false>(payload->currentMultipole,
																payload->currentEdgeLength,
																payload->currentLevel);
			break;
        } /* M2LFourier */
        case Dummy: {
            // do nothing, only serves for synchronization
            break;
        } /* Dummy */
        default:
            global_log->error() << "Undefined Quicksched task type: " << type << std::endl;
	}
}
#endif /* QUICKSCEHD */
} /* namespace bhfmm */

