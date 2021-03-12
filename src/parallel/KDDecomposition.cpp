#include "KDDecomposition.h"

#include <cfloat>
#include <sstream>
#include <fstream>
#include <climits>
#include <cmath>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "Domain.h"
#include "KDNode.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"
#include "particleContainer/adapter/FlopCounter.h"
#include "parallel/NeighbourCommunicationScheme.h"
#include "parallel/HaloRegion.h"
#include "WrapOpenMP.h"
#include "plugins/VectorizationTuner.h"

#include "KDDStaticValues.h"
#include "KDNode.h"
#include "ParticleData.h"

using namespace std;
using Log::global_log;

KDDecomposition::KDDecomposition(double cutoffRadius, int numParticleTypes, int updateFrequency,
								 int fullSearchThreshold)
	: _frequency(updateFrequency),
	  _fullSearchThreshold(fullSearchThreshold),
	  _numParticleTypes{numParticleTypes},
	  _partitionRank{calculatePartitionRank()},
	  _vecTunParticleNums(_numParticleTypes, 50) {
	_loadCalc = new TradLoad();
	_measureLoadCalc = nullptr;

	_cutoffRadius = cutoffRadius;
}

void KDDecomposition::init(Domain* domain){
    int lowCorner[KDDIM]{};
    int highCorner[KDDIM]{};
    bool coversWholeDomain[KDDIM];
    _globalNumCells = 1;

    for (int dim = 0; dim < KDDIM; dim++) {
		_globalCellsPerDim[dim] =
			static_cast<int>(floor(domain->getGlobalLength(dim) / _cutoffRadius * _cellsInCutoffRadius));
		_globalNumCells *= _globalCellsPerDim[dim];
        highCorner[dim] = _globalCellsPerDim[dim] - 1;
        _cellSize[dim] = domain->getGlobalLength(dim) / ((double) _globalCellsPerDim[dim]);
        coversWholeDomain[dim] = true;
    }

    _numParticlesPerCell.resize(_numParticleTypes * _globalNumCells);

    // create initial decomposition
    // ensure that enough cells for the number of procs are available
    _decompTree = new KDNode(_numProcs, lowCorner, highCorner, 0, 0, coversWholeDomain, 0);
    if (!_decompTree->isResolvable()) {
        auto minCellCountPerProc = std::pow(KDDStaticValues::minNumCellsPerDimension, 3);
        global_log->error() << "KDDecomposition not possible. Each process needs at least " << minCellCountPerProc
                            << " cells." << endl;
        global_log->error() << "The number of Cells is only sufficient for " << _decompTree->getNumMaxProcs() << " Procs!" << endl;
        Simulation::exit(-1);
    }
    _decompTree->buildKDTree();
    _ownArea = _decompTree->findAreaForProcess(_rank);

    // initialize the mpi data type for particles once in the beginning
    KDNode::initMPIDataType();

    global_log->info() << "Created KDDecomposition with updateFrequency=" << _frequency << ", fullSearchThreshold=" << _fullSearchThreshold << endl;

#ifdef DEBUG_DECOMP
    global_log->info() << "Initial Decomposition: " << endl;
	if (_rank == 0) {
		_decompTree->printTree("", std::cout);
	}
#endif
}

KDDecomposition::~KDDecomposition() {
//	_decompTree->serialize(string("kddecomp.dat"));
	if (_rank == 0) {
		_decompTree->plotNode("kddecomp.vtu", &_processorSpeeds);
	}
	delete _decompTree;
	KDNode::shutdownMPIDataType();
	delete _loadCalc;
	delete _measureLoadCalc;
}

void KDDecomposition::readXML(XMLfileUnits& xmlconfig) {
	/* TODO: Maybe add decomposition dimensions, default auto. */
    xmlconfig.getNodeValue("minNumCellsPerDimension", KDDStaticValues::minNumCellsPerDimension);
    global_log->info() << "KDDecomposition minNumCellsPerDimension: " << KDDStaticValues::minNumCellsPerDimension
					   << endl;
	if(KDDStaticValues::minNumCellsPerDimension==0u){
		global_log->error() << "KDDecomposition minNumCellsPerDimension has to be bigger than zero!" << std::endl;
		Simulation::exit(43);
	}
	xmlconfig.getNodeValue("updateFrequency", _frequency);
	global_log->info() << "KDDecomposition update frequency: " << _frequency << endl;
	xmlconfig.getNodeValue("fullSearchThreshold", _fullSearchThreshold);
	global_log->info() << "KDDecomposition full search threshold: " << _fullSearchThreshold << endl;
	xmlconfig.getNodeValue("heterogeneousSystems", _heterogeneousSystems);
	global_log->info() << "KDDecomposition for heterogeneous computing systems (old version, not compatible with new "
						  "VecTuner version)?: "
					   << (_heterogeneousSystems ? "yes" : "no") << endl;

   	{
		std::string deviationReductionOperation;
		xmlconfig.getNodeValue("deviationReductionOperation", deviationReductionOperation);
		if (not deviationReductionOperation.empty()) {
			if (deviationReductionOperation == "sum") {
				_deviationReductionOperation = MPI_SUM;
			} else if (deviationReductionOperation == "max") {
				_deviationReductionOperation = MPI_MAX;
			} else {
				global_log->fatal() << "Wrong deviationReductionOperation given: " << _deviationReductionOperation
									<< ". Should be 'max' or 'sum'." << std::endl;
				Simulation::exit(45681);
			}
		}
		global_log->info() << "KDDecomposition uses " << deviationReductionOperation
					   << " to reduce the deviation within the decompose step." << endl;
	}

	bool useVecTuner = false;
	xmlconfig.getNodeValue("useVectorizationTuner", useVecTuner);
	global_log->info() << "KDDecomposition using vectorization tuner: " << (useVecTuner?"yes":"no") << endl;
	if (useVecTuner){
		delete _loadCalc;
		_loadCalc = new TunerLoad();
	}

	xmlconfig.getNodeValue("clusterHetSys", _clusteredHeterogeneouseSystems);
	global_log->info() << "KDDecomposition for clustered heterogeneous systems?: " << (_clusteredHeterogeneouseSystems?"yes":"no") << endl;
	//TODO remove this check if the heterogenous Decomposition is updated to the vectorization tuner.
	if(_heterogeneousSystems){
		global_log->warning() << "The old version of the heterogeneous KDDecomposition shouldn't be used with the vectorization tuner!" << endl;
	}
	xmlconfig.getNodeValue("splitBiggestDimension", _splitBiggest);
	global_log->info() << "KDDecomposition splits along biggest domain?: " << (_splitBiggest?"yes":"no") << endl;
	xmlconfig.getNodeValue("forceRatio", _forceRatio);
	global_log->info() << "KDDecomposition forces load/performance ratio?: " << (_forceRatio?"yes":"no") << endl;

	xmlconfig.getNodeValue("rebalanceLimit", _rebalanceLimit);
	if(_rebalanceLimit > 0) {
		global_log->info() << "KDDecomposition automatic rebalancing: enabled" << endl;
		global_log->info() << "KDDecomposition rebalance limit: " << _rebalanceLimit << endl;
	}
	else {
		global_log->info() << "KDDecomposition automatic rebalancing: disabled" << endl;
	}
	xmlconfig.getNodeValue("splitThreshold", _splitThreshold);
	if(!_splitBiggest){
		global_log->info() << "KDDecomposition threshold for splitting not only the biggest Domain: " << _splitThreshold << endl;
	}

	/*
	 * Reads the Vectorization tuner parameters
	 */
	for(int i = 0; i < _numParticleTypes; ++i){
		xmlconfig.getNodeValue("particleCount" + std::to_string(i+1), _vecTunParticleNums.at(i));
		global_log->info() << "Maximum particle count in the vectorization tuner of type " << i+1 << ": " << _vecTunParticleNums.at(i) << endl;
	}
	xmlconfig.getNodeValue("generateNewFiles", _generateNewFiles);
	global_log->info() << "Generate new vectorization tuner files: " << (_generateNewFiles?"yes":"no") << endl;
	xmlconfig.getNodeValue("useExistingFiles", _useExistingFiles);
	global_log->info() << "Use existing vectorization tuner files (if available)?: " << (_useExistingFiles?"yes":"no") << endl;
	xmlconfig.getNodeValue("vecTunerAllowMPIReduce", _vecTunerAllowMPIReduce);
	global_log->info() << "Allow an MPI Reduce for the vectorization tuner?: " << (_vecTunerAllowMPIReduce?"yes":"no") << endl;

	xmlconfig.getNodeValue("doMeasureLoadCalc", _doMeasureLoadCalc);
	global_log->info() << "Use measureLoadCalc? (requires compilation with armadillo): " << (_doMeasureLoadCalc?"yes":"no") << endl;

	xmlconfig.getNodeValue("measureLoadInterpolationStartsAt", _measureLoadInterpolationStartsAt);
	global_log->info() << "measureLoad: interpolation starts at "
	                   << _measureLoadInterpolationStartsAt << endl;

	xmlconfig.getNodeValue("measureLoadIncreasingTimeValues", _measureLoadIncreasingTimeValues);
	global_log->info() << "measureLoad: Ensure that cells with more particles take longer ? "
					   << (_measureLoadIncreasingTimeValues ? "yes" : "no") << endl;

	DomainDecompMPIBase::readXML(xmlconfig);

	string oldPath(xmlconfig.getcurrentnodepath());

	xmlconfig.changecurrentnode("../datastructure");
	xmlconfig.getNodeValue("cellsInCutoffRadius", _cellsInCutoffRadius);
	global_log->info() << "KDDecomposition using cellsInCutoffRadius: " << _cellsInCutoffRadius << endl;

	// reset path
	xmlconfig.changecurrentnode(oldPath);
}

void KDDecomposition::prepareNonBlockingStage(bool /*forceRebalancing*/,
		ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber) {
	const bool removeRecvDuplicates = true;
	if(sendLeavingWithCopies()){
		DomainDecompMPIBase::prepareNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_AND_HALO_COPIES, removeRecvDuplicates);
	} else {
		DomainDecompMPIBase::prepareNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_ONLY, removeRecvDuplicates);
		DomainDecompMPIBase::prepareNonBlockingStageImpl(moleculeContainer, domain, stageNumber, HALO_COPIES, removeRecvDuplicates);
	}
}

void KDDecomposition::finishNonBlockingStage(bool /*forceRebalancing*/,
		ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber) {
	const bool removeRecvDuplicates = true;
	if(sendLeavingWithCopies()){
		DomainDecompMPIBase::finishNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_AND_HALO_COPIES, removeRecvDuplicates);
	} else {
		DomainDecompMPIBase::finishNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_ONLY, removeRecvDuplicates);
		DomainDecompMPIBase::finishNonBlockingStageImpl(moleculeContainer, domain, stageNumber, HALO_COPIES, removeRecvDuplicates);
	}
}

//check whether or not to do rebalancing in the specified step
bool doRebalancing(bool forceRebalancing, bool needsRebalance, size_t steps, int frequency){
	return forceRebalancing or ((steps % frequency == 0 or steps <= 1) and needsRebalance);
}

bool KDDecomposition::queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* /*moleculeContainer*/, Domain* /*domain*/, double etime){
	bool needsRebalance = checkNeedRebalance(etime);
	return not doRebalancing(forceRebalancing, needsRebalance, _steps, _frequency);
}

bool KDDecomposition::checkNeedRebalance(double lastTraversalTime) const {
	bool needsRebalance = false;
	if (_rebalanceLimit > 0) {
		/* automatic rebalancing */
		double localTraversalTimes[2];
		localTraversalTimes[0] = -lastTraversalTime;
		localTraversalTimes[1] = lastTraversalTime;
		double globalTraversalTimes[2];
		MPI_CHECK(MPI_Allreduce(localTraversalTimes, globalTraversalTimes, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD));
		globalTraversalTimes[0] *= -1.0;
		double timerCoeff = globalTraversalTimes[1] / globalTraversalTimes[0];
		global_log->info() << "KDDecomposition timerCoeff: " << timerCoeff << endl;
		if (timerCoeff > _rebalanceLimit) {
			needsRebalance = true;
		}
	} else {
		needsRebalance = true;
	}
	return needsRebalance;
}

void KDDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain) {
	bool needsRebalance = checkNeedRebalance(lastTraversalTime);
	bool rebalance = doRebalancing(forceRebalancing, needsRebalance, _steps, _frequency);
	_steps++;
	const bool removeRecvDuplicates = true;

	size_t measureLoadInitTimers = 2;
	if (_steps == measureLoadInitTimers and _doMeasureLoadCalc) {
		if(global_simulation->getEnsemble()->getComponents()->size() > 1){
			global_log->warning() << "MeasureLoad is designed to work with one component. Using it with more than one "
									 "component might produce bad results if their force calculation differs."
								  << std::endl;
		}
		_measureLoadCalc = new MeasureLoad(_measureLoadIncreasingTimeValues, _measureLoadInterpolationStartsAt);
	}
	size_t measureLoadStart = 50;
	if (_steps == measureLoadStart and _doMeasureLoadCalc) {
		bool faulty = _measureLoadCalc->prepareLoads(this, _comm);
		if (faulty) {
			global_log->info() << "Not using MeasureLoad as it failed. No rebalance forced." << std::endl;
		} else {
			global_log->info() << "Start using MeasureLoad, will force rebalance." << std::endl;
			delete _loadCalc;
			_loadCalc = _measureLoadCalc;
			_measureLoadCalc = nullptr;
			rebalance = true;
		}
	}

	if (not rebalance) {
		if (not moleculeContainer->isInvalidParticleReturner() or moleculeContainer->hasInvalidParticles()) {
			if (sendLeavingWithCopies()) {
				global_log->debug() << "kDD: Sending Leaving and Halos together." << std::endl;
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES,
				                                          true /*doHaloPositionCheck*/, removeRecvDuplicates);
			} else {
				global_log->debug() << "kDD: Sending Leaving, then Halos." << std::endl;
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY,
				                                          true /*doHaloPositionCheck*/, removeRecvDuplicates);
#ifndef MARDYN_AUTOPAS
				moleculeContainer->deleteOuterParticles();
#endif
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES,
				                                          true /*doHaloPositionCheck*/, removeRecvDuplicates);
			}
		} else {
			global_log->debug() << "kDD: Sending Halos." << std::endl;
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES, false /*dohaloPositionCheck*/);
		}
	} else {
		global_log->info() << "KDDecomposition: rebalancing..." << endl;
		if(moleculeContainer->isInvalidParticleReturner() and not moleculeContainer->hasInvalidParticles()){
			moleculeContainer->forcedUpdate();
		}
		if (_steps != 1) {
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY,
													  true /*doHaloPositionCheck*/, removeRecvDuplicates);
		}
		moleculeContainer->deleteOuterParticles();

		KDNode * newDecompRoot = nullptr;
		KDNode * newOwnLeaf = nullptr;

		calcNumParticlesPerCell(moleculeContainer);
		constructNewTree(newDecompRoot, newOwnLeaf, moleculeContainer);
		bool migrationSuccessful = migrateParticles(*newDecompRoot, *newOwnLeaf, moleculeContainer, domain);
		if (not migrationSuccessful) {
			global_log->error() << "A problem occurred during particle migration between old decomposition and new decomposition of the KDDecomposition." << endl;
			global_log->error() << "Aborting. Please save your input files and last available checkpoint and contact TUM SCCS." << endl;
			Simulation::exit(1);
		}
		delete _decompTree;
		_decompTree = newDecompRoot;
//		delete _ownArea; dont delete! this is a pointer only to one of the objects in the whole tree, not a real object
		_ownArea = newOwnLeaf;
		initCommunicationPartners(_cutoffRadius, domain, moleculeContainer);

		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES, true /*doHaloPositionCheck*/, removeRecvDuplicates);
	}
}

void KDDecomposition::initCommunicationPartners(double cutoffRadius, Domain * domain, ParticleContainer* moleculeContainer) {
	_neighbourCommunicationScheme->initCommunicationPartners(cutoffRadius, domain, this, moleculeContainer);
}

void KDDecomposition::getCellBorderFromIntCoords(double * lC, double * hC, int lo[3], int hi[3]) const {
	for(int d = 0 ; d < 3; ++d) {
		lC[d] = lo[d] * _cellSize[d];
		hC[d] = (hi[d] + 1) * _cellSize[d];
	}
}

void KDDecomposition::getCellIntCoordsFromRegionPeriodic(int* lo, int* hi, const double lC[3], const double hC[3], const Domain* domain) const {
	for (int d = 0; d < 3; ++d) {
		lo[d] = round(lC[d] / _cellSize[d]);
		hi[d] = round(hC[d] / _cellSize[d]) - 1; // no -1
	}
}

bool KDDecomposition::migrateParticles(const KDNode& newRoot, const KDNode& newOwnLeaf, ParticleContainer* moleculeContainer, Domain* domain) {
	// 1. compute which processes we will receive from
	// 2. issue Irecv calls
	// 3. compute which processes we will send to
	// 4. issue Isend calls
	// 5. get all

	vector<CommunicationPartner> recvPartners;
	recvPartners.clear();
	int numProcsRecv;

	vector<Molecule*> migrateToSelf;
	bool willMigrateToSelf = false;

	// issue Recv calls
	{
		// process-leaving particles have been handled, so we only need actual area
		vector<int> ranks;
		vector<int> indices;
		_decompTree->getOwningProcs(newOwnLeaf._lowCorner, newOwnLeaf._highCorner, ranks, indices);

		vector<int> numMolsToRecv;
		auto indexIt = indices.begin();
		numProcsRecv = ranks.size(); // value may change from ranks.size(), see "numProcsSend--" below
		recvPartners.reserve(numProcsRecv);
		for (unsigned i = 0; i < ranks.size(); ++i) {

			int low[3];
			int high[3];
			for (int d=0; d < 3; ++d) {
				low[d] = *(indexIt++);
				high[d] = *(indexIt++);
			}
			int numMols = 0;
			for (int iz = low[2]; iz <= high[2]; ++iz) {
				for (int iy = low[1]; iy <= high[1]; ++iy) {
					for (int ix = low[0]; ix <= high[0]; ++ ix) {
						numMols += _numParticlesPerCell[(iz * _globalCellsPerDim[1] + iy) * _globalCellsPerDim[0] + ix];
						if(_numParticleTypes >= 2){
							numMols += _numParticlesPerCell[_globalNumCells + (iz * _globalCellsPerDim[1] + iy) * _globalCellsPerDim[0] + ix];
						}
					}
				}
			}

			int partnerRank = ranks[i];

			if (partnerRank != _rank) {
				recvPartners.emplace_back(partnerRank);
				recvPartners.back().initRecv(numMols, _comm, _mpiParticleType);
			} else {
				migrateToSelf.reserve(numMols);
				// decrement numProcsRecv for following uses
				willMigrateToSelf = true;
				numProcsRecv--;
			}
		}
	}

	vector<CommunicationPartner> sendPartners;
	sendPartners.clear();
	int numProcsSend;
	// issue Send calls
	{
		// process-leaving particles have been handled, so we only need actual area
		vector<int> ranks;
		vector<int> indices;
		newRoot.getOwningProcs(_ownArea->_lowCorner, _ownArea->_highCorner, ranks, indices);

		auto indexIt = indices.begin();
		numProcsSend = ranks.size(); // value may change from ranks.size(), see "numProcsSend--" below
		sendPartners.reserve(numProcsSend);
		std::vector<Molecule> dummy;
		for (unsigned i = 0; i < ranks.size(); ++i) {
			int low[3];
			int high[3];
			double leavingLow[3];
			double leavingHigh[3];
			for (int d=0; d < 3; ++d) {
				low[d] = *(indexIt++);
				high[d] = *(indexIt++);
				leavingLow[d] = low[d] * _cellSize[d];
				leavingHigh[d] = (high[d] + 1) * _cellSize[d];
			}
			int partnerRank = ranks[i];
			if (partnerRank != _rank) {
				sendPartners.emplace_back(partnerRank, leavingLow, leavingHigh);
				const bool removeFromContainer = true;
				sendPartners.back().initSend(moleculeContainer, _comm, _mpiParticleType, LEAVING_ONLY, dummy,
											 /*don't use invalid particles*/ false, true /*do halo position change*/,
											 removeFromContainer);
				// molecules are taken out of container
			} else {
				bool inHaloRegion = true;
				for (unsigned int dimindex = 0; dimindex <3; dimindex ++){
					inHaloRegion &= leavingLow[dimindex] < getBoundingBoxMax(dimindex, domain);
					inHaloRegion &= leavingHigh[dimindex] >= getBoundingBoxMin(dimindex, domain);
				}
				if (inHaloRegion) {
					collectMoleculesInRegion(moleculeContainer, leavingLow, leavingHigh, migrateToSelf);
				}

				// decrement numProcsSend for further uses:
				mardyn_assert(willMigrateToSelf);
				numProcsSend--;
			}
		}
	}
	mardyn_assert(moleculeContainer->getNumberOfParticles() == 0ul);
	double newBoxMin[3];
	double newBoxMax[3];
	for (int dim = 0; dim < 3; dim++) {
		newBoxMin[dim] = (newOwnLeaf._lowCorner[dim]) * _cellSize[dim];
		newBoxMax[dim] = (newOwnLeaf._highCorner[dim] + 1) * _cellSize[dim];
		// for the last process the boxmax always has to be exactly domain->getGlobalLength().
		if (newOwnLeaf._highCorner[dim] + 1 == _globalCellsPerDim[dim]) {
			newBoxMax[dim] = domain->getGlobalLength(dim);
		}
	}
	bool sendTogether = moleculeContainer->rebuild(newBoxMin, newBoxMax);

	// the indirect neighborcommunicationscheme in combination with the kddecomposition is not allowed
	// to send halo and leaving particles together, as long as halo particles are not send with all data (velocity, etc.)
	bool neighborschemeAllowsDirect =
		dynamic_cast<DirectNeighbourCommunicationScheme*>(_neighbourCommunicationScheme) != nullptr;
	sendTogether &= neighborschemeAllowsDirect;
	updateSendLeavingWithCopies(sendTogether);

	global_log->set_mpi_output_all();
	double waitCounter = 5.0;
	double deadlockTimeOut = 360.0;
	bool allDone = false;
	double startTime = MPI_Wtime();
	bool migrateToSelfDone = not willMigrateToSelf;

	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (int i = 0; i < numProcsSend; ++i) {
			allDone &= sendPartners[i].testSend();
		}

		if (not migrateToSelfDone) {
			const int numMolsMigToSelf = migrateToSelf.size();
			for (int i = 0; i < numMolsMigToSelf; i++) {
				moleculeContainer->addParticle(*migrateToSelf[i]);
				delete migrateToSelf[i];
			}
			migrateToSelfDone = true;
		}

		// unpack molecules
		for (int i = 0; i < numProcsRecv; ++i) {
			allDone &= recvPartners[i].testRecv(moleculeContainer, false);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			global_log->warning() << "KDDecomposition::migrateParticles: Deadlock warning: Rank " << _rank
					<< " is waiting for more than " << waitCounter << " seconds"
					<< std::endl;
			waitCounter += 1.0;
			for (int i = 0; i < numProcsSend; ++i) {
				sendPartners[i].deadlockDiagnosticSend();
			}
			for (int i = 0; i < numProcsRecv; ++i) {
				recvPartners[i].deadlockDiagnosticRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			global_log->error() << "KDDecomposition::migrateParticles: Deadlock error: Rank " << _rank
					<< " is waiting for more than " << deadlockTimeOut
					<< " seconds" << std::endl;
			for (int i = 0; i < numProcsSend; ++i) {
				sendPartners[i].deadlockDiagnosticSend();
			}
			for (int i = 0; i < numProcsRecv; ++i) {
				recvPartners[i].deadlockDiagnosticRecv();
			}
			break;
		}

	} // while not allDone

	global_log->set_mpi_output_root(0);

	moleculeContainer->update();

	int isOK = allDone;

	MPI_Allreduce(MPI_IN_PLACE, &isOK, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	if (isOK == _numProcs) {
		return true;
	} else {
		global_log->error() << "writing checkpoint to kddecomperror.restart.dat" << std::endl;
		global_simulation->getDomain()->writeCheckpoint("kddecomperror.restart.dat", moleculeContainer, this, global_simulation->getSimulationTime());
		return false;
	}
}

void KDDecomposition::fillTimeVecs(CellProcessor **cellProc){
	if(cellProc == nullptr){
		global_log->error() << "The cellProcessor was not yet set! Please reorder fillTimeVecs, so that there won't be a problem!";
		Simulation::exit(1);
	}
	auto _tunerLoadCalc = dynamic_cast<TunerLoad*>(_loadCalc);
	if(_tunerLoadCalc){
		VectorizationTuner tuner;
		tuner.init(global_simulation->getMoleculeContainer(), &global_simulation->domainDecomposition(), global_simulation->getDomain());
		mardyn_assert(cellProc && (*cellProc));
		tuner.tune(*(_simulation.getEnsemble()->getComponents()), *_tunerLoadCalc, _vecTunParticleNums, _generateNewFiles, _useExistingFiles, _vecTunerAllowMPIReduce);
	}


}

void KDDecomposition::constructNewTree(KDNode *& newRoot, KDNode *& newOwnLeaf, ParticleContainer* moleculeContainer) {

	newRoot = new KDNode(_numProcs, &(_decompTree->_lowCorner[0]), &(_decompTree->_highCorner[0]), 0, 0, _decompTree->_coversWholeDomain, 0);
	KDNode * toCleanUp = newRoot;

	updateMeanProcessorSpeeds(_processorSpeeds,_accumulatedProcessorSpeeds, moleculeContainer);
	bool result;
	if(_clusteredHeterogeneouseSystems){
		result = heteroDecompose(newRoot, newOwnLeaf, MPI_COMM_WORLD);
	} else {
		result = decompose(newRoot, newOwnLeaf, MPI_COMM_WORLD);
	}
	if (result) {
		global_log->warning() << "Domain too small to achieve a perfect load balancing" << endl;
	}

	completeTreeInfo(newRoot, newOwnLeaf);
	delete toCleanUp;
	for (int d = 0; d < 3; ++d) {
		_neighbourCommunicationScheme->setCoverWholeDomain(d, newOwnLeaf->_coversWholeDomain[d]);
	}

	global_log->info() << "KDDecomposition: rebalancing finished" << endl;

#ifndef NDEBUG
	if (_rank == 0) {
		stringstream fname;
		fname << "kddecomp_" << _steps - 1 << ".vtu";
		newRoot->plotNode(fname.str(), &_processorSpeeds);
	}
#endif /* NDEBUG */

#ifdef DEBUG_DECOMP
	if (_rank == 0) {
		newRoot->printTree("", std::cout);
	}
#endif
}

void KDDecomposition::updateMeanProcessorSpeeds(std::vector<double>& processorSpeeds,
		std::vector<double>& accumulatedProcessorSpeeds, ParticleContainer* moleculeContainer) {
	// update the processor speed exactly twice (first update at preprocessor stage (no speeds known yet)
	// second update after first real measurements
	if(_processorSpeedUpdateCount > 1){
		return;
	}
	_processorSpeedUpdateCount++;

	if (_heterogeneousSystems) {
		FlopCounter fl(global_simulation->getcutoffRadius(), global_simulation->getLJCutoff());
		moleculeContainer->traverseCells(fl);
		double flopCount = fl.getMyFlopCount();
		double flopRate = flopCount / global_simulation->getAndResetOneLoopCompTime();
		if (flopRate == 0.) {  // for unit_tests and first simulation
			flopRate = 1.;
		}
		collCommInit(_numProcs);
		for (int i = 0; i < _rank; i++) {
			collCommAppendDouble(0.);
		}
		collCommAppendDouble(flopRate);
		for (int i = _rank + 1; i < _numProcs; i++) {
			collCommAppendDouble(0.);
		}
		collCommAllreduceSum();
		if (processorSpeeds.empty()) {
			processorSpeeds.resize(_numProcs);
			accumulatedProcessorSpeeds.clear();
		}
		for (int i = 0; i < _numProcs; i++) {
			processorSpeeds[i] = collCommGetDouble();
		}
		collCommFinalize();
	}
	else{
		if (processorSpeeds.empty()) {
			processorSpeeds.resize(_numProcs, 1.);
			_totalProcessorSpeed = _numProcs;
			_totalMeanProcessorSpeed = 1.;
			accumulatedProcessorSpeeds.clear();
		}
	}




	if (processorSpeeds.size() + 1 != accumulatedProcessorSpeeds.size()) {
		accumulatedProcessorSpeeds.resize(processorSpeeds.size() + 1);
	}
	accumulatedProcessorSpeeds[0] = 0.;
	for (size_t i = 0; i < processorSpeeds.size(); i++) {
		mardyn_assert(processorSpeeds[i] > 0);
		accumulatedProcessorSpeeds[i + 1] = accumulatedProcessorSpeeds[i] + processorSpeeds[i];
	}
	_totalProcessorSpeed =
			accumulatedProcessorSpeeds[processorSpeeds.size()];
	_totalMeanProcessorSpeed = _totalProcessorSpeed
			/ processorSpeeds.size();
}


//TODO: rewrite with min/max
double KDDecomposition::getBoundingBoxMin(int dimension, Domain* domain) {
	double globalLength = domain->getGlobalLength(dimension);
	double pos = (_ownArea->_lowCorner[dimension]) * _cellSize[dimension];
	if (pos < 0) {
		return 0;
	}
	else if (pos > globalLength) {
		return globalLength;
	}
	else {
		return pos;
	}
}

double KDDecomposition::getBoundingBoxMax(int dimension, Domain* domain) {
	if (_ownArea->_highCorner[dimension] + 1 == _globalCellsPerDim[dimension]) {
		return domain->getGlobalLength(dimension);
	}
	double globalLength = domain->getGlobalLength(dimension);
	double pos = (_ownArea->_highCorner[dimension] + 1) * _cellSize[dimension];
	if (pos < 0) {
		return 0;
	}
	else if (pos > globalLength) {
		return globalLength;
	}
	else {
		return pos;
	}
}

void KDDecomposition::printDecomp(const std::string& filename, Domain* domain) {
	if (_rank == 0) {
		ofstream povcfgstrm(filename.c_str());
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " " << domain->getGlobalLength(2) << endl;
		povcfgstrm << "decompData Regions" << endl;
		povcfgstrm.close();
	}

	stringstream output;
	output  << getBoundingBoxMin(0,domain) << " " << getBoundingBoxMin(1,domain) << " "
			<< getBoundingBoxMin(2,domain) << " " << getBoundingBoxMax(0,domain) << " "
			<< getBoundingBoxMax(1,domain) << " " << getBoundingBoxMax(2,domain) << "\n";
	string output_str = output.str();
#ifdef ENABLE_MPI
	MPI_File fh;
	MPI_File_open(_comm, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_APPEND | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	uint64_t write_size = output_str.size();
	uint64_t offset = 0;
	if(_rank == 0) {
		MPI_Offset file_end_pos;
		MPI_File_seek(fh, 0, MPI_SEEK_END);
		MPI_File_get_position(fh, &file_end_pos);
		write_size += file_end_pos;
		MPI_Exscan(&write_size, &offset, 1, MPI_UINT64_T, MPI_SUM, _comm);
		offset += file_end_pos;
	} else {
		MPI_Exscan(&write_size, &offset, 1, MPI_UINT64_T, MPI_SUM, _comm);
	}
	MPI_File_write_at(fh, offset, output_str.c_str(), output_str.size(), MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
#else
	ofstream povcfgstrm(filename.c_str(), ios::app);
	povcfgstrm << output_str;
	povcfgstrm.close();
#endif
}


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$ private Methoden, die von exchangeMolecule benötigt werden $
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

void KDDecomposition::completeTreeInfo(KDNode*& root, KDNode*& ownArea) {

	int numElementsToRecv = root->_numProcs * 2 - 1;
	vector<KDNode*> ptrToAllNodes;
	vector<int> child1;
	vector<int> child2;
	ptrToAllNodes.resize(numElementsToRecv);
	child1.resize(numElementsToRecv);
	child2.resize(numElementsToRecv);

	KDNode* oldNode = root;
	// each process walks down the tree until an owned element is reached
	while (oldNode->_owningProc != _rank) {
		if (_rank < oldNode->_owningProc + oldNode->_child1->_numProcs) {
			oldNode = oldNode->_child1;
		}
		else {
			oldNode = oldNode->_child2;
		}
	}

	int nextSendingProcess = 0;
	for (int nodeID = 0; nodeID < numElementsToRecv; nodeID++) {

		KDNode::MPIKDNode mpiKDNode;
		if (oldNode->_nodeID == nodeID) {
			mpiKDNode = oldNode->getMPIKDNode();

			if (oldNode->_numProcs > 1) {
				oldNode = oldNode->_child1;
			}
		}

		MPI_CHECK( MPI_Bcast(&mpiKDNode, 1, KDNode::MPIKDNode::Datatype, nextSendingProcess, MPI_COMM_WORLD) );
		bool coversAll[3];

		coversAll[0] = mpiKDNode.getCoversWholeDomain(0);
		coversAll[1] = mpiKDNode.getCoversWholeDomain(1);
		coversAll[2] = mpiKDNode.getCoversWholeDomain(2);

		ptrToAllNodes[nodeID] = new KDNode(mpiKDNode.getNumProcs(), mpiKDNode.getLowCorner(),
				mpiKDNode.getHighCorner(), mpiKDNode.getNodeID(), mpiKDNode.getOwningProc(), coversAll, mpiKDNode.getLevel());
		ptrToAllNodes[nodeID]->_load = mpiKDNode.getLoad();
		ptrToAllNodes[nodeID]->_optimalLoadPerProcess = mpiKDNode.getOptimalLoadPerProcess();
		ptrToAllNodes[nodeID]->_deviationLowerBound = mpiKDNode.getDeviationLowerBound();
		ptrToAllNodes[nodeID]->_deviation = mpiKDNode.getDeviation();
		child1[nodeID] = mpiKDNode.getFirstChildID();
		child2[nodeID] = mpiKDNode.getSecondChildID();
		nextSendingProcess = mpiKDNode.getNextSendingProcess();
	}

	// connect all Nodes of the tree
	for (int nodeID = 0; nodeID < numElementsToRecv; nodeID++) {
		// if node is not leaf node, connect to the two children
		if (child1[nodeID] >= 0) {
			ptrToAllNodes[nodeID]->_child1 = ptrToAllNodes[child1[nodeID]];
			ptrToAllNodes[nodeID]->_child2 = ptrToAllNodes[child2[nodeID]];
		}
	}

	root = ptrToAllNodes[0];
	ownArea = ptrToAllNodes[ownArea->_nodeID];
}

#ifdef DEBUG_DECOMP
void printChildrenInfo(std::ofstream& filestream, KDNode* node, double minDev) {
	for (int i = 0; i < node->_level; i++) { filestream << "   ";}
	filestream << " * " << "load=" << node->_load << " optLoad=" << node->_optimalLoadPerProcess << " expDev=" << node->_deviationLowerBound << " minDev=" << minDev << endl;
	for (int i = 0; i < node->_level; i++) { filestream << "   ";}
	filestream << "   [" << node->_child1->_lowCorner[0] << "," << node->_child1->_lowCorner[1] << "," << node->_child1->_lowCorner[2] << "]"
			<< "[" << node->_child1->_highCorner[0] << "," << node->_child1->_highCorner[1] << "," << node->_child1->_highCorner[2] << "]"
			<< " load=" << node->_child1->_load << " #procs=" << node->_child1->_numProcs << " avgLoad=" << node->_child1->calculateAvgLoadPerProc() << endl;
	for (int i = 0; i < node->_level; i++) { filestream << "   ";}
	filestream << "   [" << node->_child2->_lowCorner[0] << "," << node->_child2->_lowCorner[1] << "," << node->_child2->_lowCorner[2] << "]"
				<< "[" << node->_child2->_highCorner[0] << "," << node->_child2->_highCorner[1] << "," << node->_child2->_highCorner[2] << "]"
				<< " load=" << node->_child2->_load << " #procs=" << node->_child2->_numProcs << " avgLoad=" << node->_child2->calculateAvgLoadPerProc() << endl;
}
#endif

bool KDDecomposition::decompose(KDNode* fatherNode, KDNode*& ownArea, MPI_Comm commGroup) {
	return decompose(fatherNode, ownArea, commGroup, FLT_MAX);
}

bool KDDecomposition::decompose(KDNode* fatherNode, KDNode*& ownArea, MPI_Comm commGroup, const double globalMinimalDeviation) {
	bool domainTooSmall = false;
	// recursion termination criterion
	if (fatherNode->_numProcs == 1) {
		// own area must belong to this process!
		mardyn_assert(fatherNode->_owningProc == _rank);
		ownArea = fatherNode;
		fatherNode->calculateDeviation(&_processorSpeeds, _totalMeanProcessorSpeed);
		return domainTooSmall;
	}

	std::list<KDNode*> possibleSubdivisions;
	domainTooSmall = calculateAllPossibleSubdivisions(fatherNode, possibleSubdivisions, commGroup);
	mardyn_assert(not possibleSubdivisions.empty());

	KDNode* bestSubdivision = nullptr;
	double minimalDeviation = globalMinimalDeviation;
	auto iter = possibleSubdivisions.begin();
	int iterations = 0;
	int log2proc = 0;// calculates the logarithm of _numProcs (base 2)
	while ((_numProcs >> log2proc) > 1) {
		log2proc++;
	}
	// if we are near the root of the tree, we just take the first best subdivision
	// this is even valid for heterogeneous balancing, since the possibleSubdivisions are calculated accordingly, i.e., they are
	// sorted by their expected deviation (lowest to biggest).
	int maxIterations = 1;
	if (fatherNode->_level > (log2proc - _fullSearchThreshold)) {//only do proper search, if this condition is fulfilled (numProcs < 2^(level + threshold))
		maxIterations = INT_MAX;
	}

#ifdef DEBUG_DECOMP
	std::stringstream fname;
	fname << "Div_proc_" << _rank << "_step_" << _steps << ".txt";
	std::ofstream filestream(fname.str().c_str(), ios::app);
	filestream.precision(8);
	for (int i = 0; i < fatherNode->_level; i++) { filestream << "   ";}
	filestream << "Division at rank=" << _rank << " for [" << fatherNode->_lowCorner[0]
               << ","<< fatherNode->_lowCorner[1] << "," << fatherNode->_lowCorner[2] << "] [" << fatherNode->_highCorner[0]
               << ","<< fatherNode->_highCorner[1] << "," << fatherNode->_highCorner[2] << "] " <<
               "level=" << fatherNode->_level << " #divisions=" << possibleSubdivisions.size() << endl;
#endif

	while (iter != possibleSubdivisions.end() && (iterations < maxIterations) && (*iter)->_deviationLowerBound < minimalDeviation) {
		iterations++;
#ifdef DEBUG_DECOMP
		printChildrenInfo(filestream, *iter, minimalDeviation);
#endif

		// compute the next subdivision depending on the current rank (either first or second subdivision)
		vector<int> origRanks;
		int newNumProcs;
		if (_rank < (*iter)->_child2->_owningProc) {  // assign the current rank to either the first (child1)...
			origRanks.resize((*iter)->_child1->_numProcs);
			for (int i = 0; i < (*iter)->_child1->_numProcs; i++) {
				origRanks[i] = i;  // this group will consist of the first (*iter)->_child1->_numProcs processes/ranks of the current communicator (origGroup)
			}
			newNumProcs = (*iter)->_child1->_numProcs;
		}
		else {                                       // ...or the second group (child2) for the MPI communication
			origRanks.resize((*iter)->_child2->_numProcs);
			for (int i = 0; i < (*iter)->_child2->_numProcs; i++) {
				origRanks[i] = i + (*iter)->_child1->_numProcs; // this group consists of the last (*iter)->_child2->_numProcs processes/ranks of the current communicator (origGroup)
			}
			newNumProcs = (*iter)->_child2->_numProcs;
		}

		MPI_Comm newComm;
		MPI_Group origGroup, newGroup;

		MPI_CHECK( MPI_Comm_group(commGroup, &origGroup) );
		MPI_CHECK( MPI_Group_incl(origGroup, newNumProcs, &origRanks[0], &newGroup) );//create new MPI group based on rank (as calculated before)
		MPI_CHECK( MPI_Comm_create(commGroup, newGroup, &newComm) );

		KDNode* newOwnArea = nullptr;
		double deviationChildren[] = {0.0, 0.0};

		if (_rank < (*iter)->_child2->_owningProc) {  // compute the subdivision of the first child ...
			// do not use the function call directly in the logical expression, as it may
			// not be executed due to conditional / short-circuit evaluation!
			bool subdomainTooSmall = decompose((*iter)->_child1, newOwnArea, newComm, minimalDeviation);
			deviationChildren[0] = (*iter)->_child1->_deviation;
			domainTooSmall = (domainTooSmall || subdomainTooSmall);
		} else {									  // ... or the second child
			mardyn_assert(_rank >= (*iter)->_child2->_owningProc);
			bool subdomainTooSmall = decompose((*iter)->_child2, newOwnArea, newComm, minimalDeviation);
			deviationChildren[1] = (*iter)->_child2->_deviation;
			domainTooSmall = (domainTooSmall || subdomainTooSmall);
		}

		MPI_CHECK( MPI_Group_free(&newGroup));
		MPI_CHECK( MPI_Comm_free(&newComm) );
		MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, deviationChildren, 2, MPI_DOUBLE, _deviationReductionOperation, commGroup));  // reduce the deviations
		(*iter)->_child1->_deviation = deviationChildren[0];
		(*iter)->_child2->_deviation = deviationChildren[1];
		(*iter)->calculateDeviation();
		if((*iter)->_deviation < (*iter)->_deviationLowerBound){
			global_log->warning() << "Calculated deviation " << (*iter)->_deviation << " lower than lower bound "
					  << (*iter)->_deviationLowerBound << ". This should not happen. Please report a bug." << std::endl;
		}
#ifdef DEBUG_DECOMP
		for (int i = 0; i < fatherNode->_level; i++) { filestream << "   ";}
		filestream << "   deviation=" << (*iter)->_deviation << " (ch1:" << deviationChildren[0] << "ch2:" << deviationChildren[1] << endl;
#endif
		if ((*iter)->_deviation < minimalDeviation) {
			delete bestSubdivision;// (deleting of nullptr is ok and does not produce errors, since delete checks for nullptr)
			bestSubdivision = *iter;
			minimalDeviation = (*iter)->_deviation;
			ownArea = newOwnArea;
		} else {
			delete *iter;
		}
		iter++;
	}

	while (iter != possibleSubdivisions.end()) {
		delete *iter;
		iter++;
	}

	// reassign children and delete cloned node, if a solution
	// was found in this subtree.
	if (bestSubdivision == nullptr) {
		fatherNode->_deviation = FLT_MAX;
	} else {
		*fatherNode = *bestSubdivision;  // assignment operator (NOT copy operator) -> also assigns children to fatherNode
		bestSubdivision->_child1 = nullptr;  // remove children from bestSubdivision, otherwise they will be deleted
		bestSubdivision->_child2 = nullptr;  // remove children from bestSubdivision, otherwise they will be deleted
		delete bestSubdivision;
	}
#ifdef DEBUG_DECOMP
	filestream.close();
#endif
	return domainTooSmall;
}


bool KDDecomposition::calculateAllPossibleSubdivisions(KDNode* node, std::list<KDNode*>& subdividedNodes, MPI_Comm commGroup) {
	bool domainTooSmall = false;
	vector<vector<double> > costsLeft(3);
	vector<vector<double> > costsRight(3);
	calculateCostsPar(node, costsLeft, costsRight, commGroup);
	global_log->debug() << "calculateAllPossibleSubdivisions: " << std::endl;
	double leftRightLoadRatio = 1.;  // divide load 50/50 -> this can be changed later on if the topology of the system should be taken into account.
	// update leftRightRatio to something meaningful (that is representable by the compute power distribution of the processes:

	// speed that is close to LeftRightRatio
	double optimalSpeed = (_accumulatedProcessorSpeeds[node->_owningProc + node->_numProcs]
			- _accumulatedProcessorSpeeds[node->_owningProc]) * leftRightLoadRatio / (1. + leftRightLoadRatio);
	double searchSpeed = optimalSpeed + _accumulatedProcessorSpeeds[node->_owningProc];
	{
	auto iter = std::lower_bound(_accumulatedProcessorSpeeds.begin() + node->_owningProc,
			_accumulatedProcessorSpeeds.begin() + node->_owningProc + node->_numProcs + 1, searchSpeed); // +1 since _accumulatedProcessorSpeeds are shifted and of size (numprocs+1)
	// returns iterator to first instance of array, that is >= optimalSpeed = totalspeed * rho / (1+rho)
	// calculated as following: rho * R = L; L leftLoad, R rightLoad, rho=ratio
	// G = L + R   => rho * (G - L) = L   => L = rho / (1 + rho) * G

	int leftRightLoadRatioIndex;
	if (fabs(*(iter - 1) - searchSpeed) < fabs(*iter - searchSpeed)) { // iter-1 always exists, since _accumulatedProcessorSpeeds[0] = 0
		leftRightLoadRatioIndex = min((int) (iter - 1 - _accumulatedProcessorSpeeds.begin() - node->_owningProc), node->_numProcs - 1);
	} else {
		leftRightLoadRatioIndex = min((int) (iter - _accumulatedProcessorSpeeds.begin() - node->_owningProc), node->_numProcs - 1);
	}
	leftRightLoadRatioIndex = max(1, leftRightLoadRatioIndex);

	leftRightLoadRatio = (_accumulatedProcessorSpeeds[node->_owningProc + leftRightLoadRatioIndex] - _accumulatedProcessorSpeeds[node->_owningProc]) /
			(_accumulatedProcessorSpeeds[node->_owningProc + node->_numProcs] - _accumulatedProcessorSpeeds[node->_owningProc + leftRightLoadRatioIndex]);

	}

	bool splitLoad = true;  // indicates, whether to split the domain according to the load
							// or whether the domain should simply be split in half and the number of processes should be distributed accordingly.

	size_t dimInit = 0;
	size_t dimEnd = 3;
	if(_splitBiggest or _splitThreshold < node->_numProcs){
		size_t max = costsLeft[0].size();
		size_t maxInd = 0;
		for (unsigned int dim = 1; dim < 3; dim++){
			if (costsLeft[dim].size() > max){
				max = costsLeft[dim].size();
				maxInd = dim;
			}
		}
		dimInit = maxInd;
		dimEnd = dimInit + 1;
	}

	for (size_t dim = dimInit; dim < dimEnd; dim++) {
		// we need at least 2*KDDStaticValues::minNumCellsPerDimension cells in this direction (= 2*KDDStaticValues::minNumCellsPerDimension different splitting planes)
		if (costsLeft[dim].size() < 2 * KDDStaticValues::minNumCellsPerDimension) {
			continue;
		}
		// if a node has some more processors, we probably don't have to find the
		// "best" partitioning, but it's sufficient to divide the node in the middle
		// and shift the processor count according to the load imbalance.

		// loop only from KDDStaticValues::minNumCellsPerDimension - 1 to max-
		// (KDDStaticValues::minNumCellsPerDimension - 1) (instead 0 to max) to avoid cell-regions with size <
		// KDDStaticValues::minNumCellsPerDimension
		int startIndex = KDDStaticValues::minNumCellsPerDimension - 1;
		int maxEndIndex =
			node->_highCorner[dim] - node->_lowCorner[dim] - (KDDStaticValues::minNumCellsPerDimension - 1);
		int endIndex = maxEndIndex;

		if (node->_numProcs > _fullSearchThreshold or _forceRatio) {
			if (splitLoad) {  // we choose the index to be the best possible for splitting the ratios.
				double minError = fabs(costsLeft[dim][0] / costsRight[dim][0] - leftRightLoadRatio);
				size_t index = 0;
				double error;
				for (size_t i = 1; i < costsLeft[dim].size(); ++i) {
					error = fabs(costsLeft[dim][i] / costsRight[dim][i] - leftRightLoadRatio);
					if (error < minError) {
						minError = error;
						index = i;
					}
				}
				startIndex = max((size_t)startIndex, index);
				endIndex = min(startIndex + 1, endIndex);

				global_log->debug() << "splitLoad: startindex " << index << " of " << costsLeft[dim].size() <<std::endl;
			} else {  // If we have more than _fullSearchThreshold processes left, we split the domain in half.
				startIndex = max(startIndex, (node->_highCorner[dim] - node->_lowCorner[dim] - 1) / 2);
				endIndex = min(endIndex, startIndex + 1);
			}
		}


		int i = startIndex;
		while ( ((i < endIndex) || (subdividedNodes.empty() && i < maxEndIndex))) {
			//for (int i = startIndex; i < endIndex; i++) {

			double optCostPerProc = (costsLeft[dim][i] + costsRight[dim][i]) / ((double) node->_numProcs);
			int optNumProcsLeft;
			if (splitLoad) {  // if we split the load in a specific ratio, numProcsLeft is calculated differently
				if(_accumulatedProcessorSpeeds.empty()){
					global_log->error() << "no processor speeds given" << std::endl;
					Simulation::exit(-1);
				}
				double optimalLoad = (_accumulatedProcessorSpeeds[node->_owningProc + node->_numProcs] - _accumulatedProcessorSpeeds[node->_owningProc]) * leftRightLoadRatio
						/ (1. + leftRightLoadRatio);
				double searchLoad = optimalLoad + _accumulatedProcessorSpeeds[node->_owningProc];
				auto iter = std::lower_bound(_accumulatedProcessorSpeeds.begin() + node->_owningProc,
						_accumulatedProcessorSpeeds.begin() + node->_owningProc + node->_numProcs + 1, searchLoad);  // +1 since _accumulatedProcessorSpeeds are shifted and of size (numprocs+1)
				// returns iterator to first instance of array, that is >= optimalSpeed = totalspeed * rho / (1+rho)
				// calculated as following: rho * R = L; L leftLoad, R rightLoad, rho=ratio
				// G = L + R   => rho * (G - L) = L   => L = rho / (1 + rho) * G

				if (fabs(*(iter - 1) - searchLoad) < fabs(*iter - searchLoad)) {  // iter-1 always exists, since _accumulatedProcessorSpeeds[0] = 0
					optNumProcsLeft = min((int)(iter - 1 - _accumulatedProcessorSpeeds.begin() - node->_owningProc), node->_numProcs - 1);
				} else {
					optNumProcsLeft = min((int)(iter - _accumulatedProcessorSpeeds.begin() - node->_owningProc), node->_numProcs - 1);
				}

			} else{
				optNumProcsLeft = min(round(costsLeft[dim][i] / optCostPerProc), (double) (node->_numProcs - 1));
			}

			int numProcsLeft = max(1, optNumProcsLeft);

			auto* clone = new KDNode(*node);
			if (clone->_level == 0) {
				clone->_optimalLoadPerProcess = optCostPerProc;
			}

			clone->split(dim, node->_lowCorner[dim] + i, numProcsLeft);
			if ( (unsigned int) (clone->_child1->_numProcs + clone->_child2->_numProcs) >
			        (clone->_child1->getNumMaxProcs() + clone->_child2->getNumMaxProcs())) {
				domainTooSmall = true;
				delete clone;
				// domain is not resolvable at all -> continue
				// I think, this case should not happen, but Martin Buchholz coded it in his code...
				// OK, it happens...!
				i++;
				continue;  // no break, since other configurations can be proper (due to division of whole numbers (ganzzahldivision))
			}

			while ( (! clone->_child1->isResolvable()) && clone->_child2->isResolvable()) {
				// shift procs to child 2, adapt owner of child 2
				clone->_child1->_numProcs--;
				clone->_child2->_numProcs++;
				clone->_child2->_owningProc--;
				clone->_child2->_nodeID = clone->_nodeID + 2 * clone->_child1->_numProcs;
				domainTooSmall = true;
			}

			while ( clone->_child1->isResolvable() && (! clone->_child2->isResolvable())) {
				// shift procs to child 1, , adapt owner of child 2
				clone->_child1->_numProcs++;
				clone->_child2->_numProcs--;
				clone->_child2->_owningProc++;
				clone->_child2->_nodeID = clone->_nodeID + 2 * clone->_child1->_numProcs;
				domainTooSmall = true;
			}

			// This should never happen, but it just means, that the domain could not be split properly.
			if ((clone->_child1->_numProcs <= 0 || clone->_child1->_numProcs >= node->_numProcs) ||
					(clone->_child2->_numProcs <= 0 || clone->_child2->_numProcs >= node->_numProcs) ){
				//continue;
				global_log->error_always_output() << "ERROR in calculateAllPossibleSubdivisions(), part of the domain was not assigned to a proc" << endl;
				Simulation::exit(1);
			}
			mardyn_assert( clone->_child1->isResolvable() && clone->_child2->isResolvable() );

			clone->_child1->_load = costsLeft[dim][i];
			clone->_child2->_load = costsRight[dim][i];
			clone->_load = costsLeft[dim][i] + costsRight[dim][i];
			clone->calculateDeviationLowerBound(&_accumulatedProcessorSpeeds);

			// sort node according to expected deviation
			auto iter = subdividedNodes.begin();
			while (iter != subdividedNodes.end() && ((*iter)->_deviationLowerBound < clone->_deviationLowerBound)) {
				iter++;
			}
			subdividedNodes.insert(iter, clone);
			i++;
		}
	}
	return domainTooSmall;
}

/**
 * Calculate the division cost for all possible divisions of the node area.
 *
 * Is done in the following way:
 * - get the number of particles per cell globally (i.e. for all cell in the domain)
 * - get the number of particle pairs per cell globally
 * - then calculate the costs for all possible subdivisions.
 */
void KDDecomposition::calculateCostsPar(KDNode* area, vector<vector<double> >& costsLeft, vector<vector<double> >& costsRight, MPI_Comm commGroup) {
	vector<vector<double> > cellCosts;
	cellCosts.resize(3);

	int newRank;
	MPI_Group group;

	MPI_CHECK( MPI_Comm_group(commGroup, &group) );
	MPI_CHECK( MPI_Group_rank(group, &newRank) );

	int dimStartIndex[3];
	int dimStopIndex[3];
	dimStartIndex[0] = 0;
	dimStopIndex[0] = (area->_highCorner[0] - area->_lowCorner[0]);
	for (int dim = 1; dim < 3; dim++) {
		dimStartIndex[dim] = dimStopIndex[dim - 1] + 1;
		dimStopIndex[dim] = dimStopIndex[dim - 1] + (area->_highCorner[dim] - area->_lowCorner[dim] + 1);
	}
	int sumNumLayers = dimStopIndex[2] + 1;
	int startIndex = sumNumLayers * newRank / area->_numProcs;
	int stopIndex = sumNumLayers * (newRank + 1) / area->_numProcs - 1;
	for (int dim = 0; dim < 3; dim++) {

		cellCosts[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		costsLeft[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		costsRight[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);

		// if this process doesn't has to do anything in this dimension, continue
		if (startIndex > dimStopIndex[dim] || stopIndex < dimStartIndex[dim]) continue;

		int loopstart = max(0, startIndex - dimStartIndex[dim]);
		int loopend = min(area->_highCorner[dim] - area->_lowCorner[dim], (area->_highCorner[dim] - area->_lowCorner[dim]) + stopIndex - dimStopIndex[dim]);

		bool sendCostValue = false;
		bool recvCostValue = false;
		if (loopstart > 0)
			recvCostValue = true;
		if (loopend < area->_highCorner[dim] - area->_lowCorner[dim])
			sendCostValue = true;

		for (int i_dim = loopstart; i_dim <= loopend; i_dim++) {

			if (i_dim == 0)
				cellCosts[dim][i_dim] = 0;
			else
				cellCosts[dim][i_dim] = cellCosts[dim][i_dim - 1];

			int dim1, dim2;
			if (dim == 0) {
				dim1 = 1;
				dim2 = 2;
			}
			else if (dim == 1) {
				dim1 = 0;
				dim2 = 2;
			}
			else {
				dim1 = 0;
				dim2 = 1;
			}
			for (int i_dim1 = 0; i_dim1 <= area->_highCorner[dim1] - area->_lowCorner[dim1]; i_dim1++) {
				for (int i_dim2 = 0; i_dim2 <= area->_highCorner[dim2] - area->_lowCorner[dim2]; i_dim2++) {

					int numParts1 = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, i_dim, i_dim1, i_dim2, area)];
					int numParts2 = _numParticleTypes == 1 ? 0 :
							_numParticlesPerCell[_globalNumCells + getGlobalIndex(dim, dim1, dim2, i_dim, i_dim1, i_dim2, area)];

					//_maxPars = max(_maxPars, numParts);
					_maxPars = max(_maxPars, numParts1);
					_maxPars2 = max(_maxPars2, numParts2);
					// #######################
					// ## Cell Costs        ##
					// #######################
					cellCosts[dim][i_dim] += _loadCalc->getOwn(numParts1, numParts2);

					// all Neighbours
					for (int neighInd_divDim = i_dim - 1; neighInd_divDim <= i_dim + 1; neighInd_divDim++) {
						int nI_dim = neighInd_divDim;
						if (nI_dim + area->_lowCorner[dim] < 0)
							nI_dim = _globalCellsPerDim[dim] - 1;
						else if (nI_dim + area->_lowCorner[dim] >= _globalCellsPerDim[dim])
							nI_dim = 0;
						//counts if the the offset in this dimension is zero
						//if the offset in zero dimensions is zero the neighbor shares a corner
						//1 => edge
						//2 => face
						//3 => the cell itself
						const int zeroCount1 = (neighInd_divDim == i_dim) ? 1 : 0;

						for (int neighInd_dim1 = i_dim1 - 1; neighInd_dim1 <= i_dim1 + 1; neighInd_dim1++) {
							// adjust index in case of periodic boundary
							int nI_dim1 = neighInd_dim1;
							if (nI_dim1 + area->_lowCorner[dim1] < 0)
								nI_dim1 = _globalCellsPerDim[dim1] - 1;
							else if (nI_dim1 + area->_lowCorner[dim1] >= _globalCellsPerDim[dim1])
								nI_dim1 = 0;
							const int zeroCount2 = (neighInd_dim1 == i_dim1) ? 1 : 0;

							for (int neighInd_dim2 = i_dim2 - 1; neighInd_dim2 <= i_dim2 + 1; neighInd_dim2++) {
								// adjust index in case of periodic boundary
								int nI_dim2 = neighInd_dim2;
								if (nI_dim2 + area->_lowCorner[dim2] < 0)
									nI_dim2 = _globalCellsPerDim[dim2] - 1;
								else if (nI_dim2 + area->_lowCorner[dim2] >= _globalCellsPerDim[dim2])
									nI_dim2 = 0;
								const int zeroCount3 = (neighInd_dim2 == i_dim2) ? 1 : 0;
								// count only forward neighbours
								/*
								if (getGlobalIndex(dim, dim1, dim2, nI_dim, nI_dim1, nI_dim2, area) > getGlobalIndex(dim, dim1, dim2, i_dim, i_dim1, i_dim2, area)) {
									int numPartsNeigh = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, nI_dim, nI_dim1, nI_dim2, area)];
									cellCosts[dim][i_dim] += 0.5 * (double) (numParts * numPartsNeigh);
								}
								*/
								const int numPartsNeigh1 = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, nI_dim, nI_dim1, nI_dim2, area)];
								int numPartsNeigh2= _numParticleTypes == 1 ? 0 :
															_numParticlesPerCell[_globalNumCells + getGlobalIndex(dim, dim1, dim2, nI_dim, nI_dim1, nI_dim2, area)];

								switch(zeroCount1+zeroCount2+zeroCount3){
								case 0: //the current neighbour shares a corner
									cellCosts[dim][i_dim] += _loadCalc->getCorner(numParts1, numParts2);
									break;
								case 1: //edge
									cellCosts[dim][i_dim] += _loadCalc->getEdge(numParts1, numParts2);
									break;
								case 2: //face
									cellCosts[dim][i_dim] += _loadCalc->getFace(numParts1, numParts2);
									break;
									//3 zeroes is the cell itself which was already counted
								}
							}
						}
					}
				}
			}
		}

		// exchange intermediate calc costs
		MPI_Status recvStat;
		double tempRecvCosts, tempSendCosts;
		tempSendCosts = 0;
		if (recvCostValue) {
			MPI_CHECK( MPI_Recv(&tempRecvCosts, 1, MPI_DOUBLE, _rank - 1, 123, MPI_COMM_WORLD, &recvStat) );
			if (sendCostValue) {
				tempSendCosts = tempRecvCosts;
			}
		}
		if (sendCostValue) {
			tempSendCosts += cellCosts[dim][loopend];
			MPI_CHECK( MPI_Send(&tempSendCosts, 1, MPI_DOUBLE, _rank + 1, 123, MPI_COMM_WORLD) );
		}
		if (recvCostValue) {
			for (int i_dim = loopstart; i_dim <= loopend; i_dim++) {
				cellCosts[dim][i_dim] += tempRecvCosts;
			}
		}
	}
	for (int dim = 0; dim < 3; dim++) {
		vector<double> cellCostsSum;
		cellCostsSum.resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);

		int size2 = cellCostsSum.size();
		MPI_CHECK( MPI_Allreduce(&cellCosts[dim][0], &cellCostsSum[0], size2, MPI_DOUBLE, MPI_SUM, commGroup) );

		for (int i_dim = 0; i_dim <= area->_highCorner[dim] - area->_lowCorner[dim]; i_dim++) {
			costsLeft[dim][i_dim] = cellCostsSum[i_dim];
			costsRight[dim][i_dim] = cellCostsSum[area->_highCorner[dim] - area->_lowCorner[dim]] - cellCostsSum[i_dim];
		}
	}
}


unsigned int KDDecomposition::getGlobalIndex(int divDim, int dim1, int dim2, int index_divDim, int index_dim1, int index_dim2, KDNode* localArea) {
	int xIndex, yIndex, zIndex;
	if (divDim == 0) {
		xIndex = index_divDim + localArea->_lowCorner[divDim];
		yIndex = index_dim1 + localArea->_lowCorner[dim1];
		zIndex = index_dim2 + localArea->_lowCorner[dim2];
	}
	else if (divDim == 1) {
		xIndex = index_dim1 + localArea->_lowCorner[dim1];
		yIndex = index_divDim + localArea->_lowCorner[divDim];
		zIndex = index_dim2 + localArea->_lowCorner[dim2];
	}
	else {
		xIndex = index_dim1 + localArea->_lowCorner[dim1];
		yIndex = index_dim2 + localArea->_lowCorner[dim2];
		zIndex = index_divDim + localArea->_lowCorner[divDim];
	}
	return _globalCellsPerDim[0] * (zIndex * _globalCellsPerDim[1] + yIndex) + xIndex;
}


//##########################################################################
//##########################################################################
//###                                                                    ###
//###                  evtl unnoetige private Methoden                   ###
//###                                                                    ###
//##########################################################################
//##########################################################################

int KDDecomposition::ownMod(int number, int modulo) const {
	int result = number % modulo;
	if (result < 0)
		result += modulo;
	return result;
}

// TODO: this method could or should be moved to KDNode.
void KDDecomposition::getOwningProcs(int low[KDDIM], int high[KDDIM], KDNode* decompTree, KDNode* testNode, vector<int>* procIDs, vector<int>* neighbHaloAreas) const {
	// For areas overlapping the domain given by decompTree, the overlapping part is
	// mapped to the corresponding area on the other side of the domain (periodic boundary)
	// The boolean variable overlap stores for each coordinate direction whether the area overlaps.
	bool overlap[KDDIM];
	bool coversWholeDomain[KDDIM];
	for (int dim = 0; dim < KDDIM; dim++) {
		// if the domain is e.g. 10 cells long, and one proc has 0-8, halo will go from -1 to 9,
		// which equals (periodic boundary) 9 to 9, which would be the same as a proc which has onle
		// this one cell. Therefore, the flag coversWholeRegion is introduced
		// @todo Do this nicer
		if (high[dim] - low[dim] >= decompTree->_highCorner[dim] - decompTree->_lowCorner[dim])
			coversWholeDomain[dim] = true;
		else
			coversWholeDomain[dim] = false;

		low[dim] = ownMod(low[dim], (decompTree->_highCorner[dim] - decompTree->_lowCorner[dim] + 1));
		high[dim] = ownMod(high[dim], (decompTree->_highCorner[dim] - decompTree->_lowCorner[dim] + 1));

		if (low[dim] > high[dim])
			overlap[dim] = true;
		else
			overlap[dim] = false;

		if (coversWholeDomain[dim])
			low[dim] = high[dim] + 1;
	}

	// First find out whether the area intersects the area given by testNode
	bool areasIntersect = true;
	for (int dim = 0; dim < KDDIM; dim++) {
		if (not (coversWholeDomain[dim])) { // TODO: check whether this can be removed || _ownArea->_coversWholeDomain[dim])){
			if (overlap[dim]) {
				if (high[dim] < testNode->_lowCorner[dim] && low[dim] > testNode->_highCorner[dim])
					areasIntersect = false;
			}
			else {
				if (high[dim] < testNode->_lowCorner[dim] || low[dim] > testNode->_highCorner[dim])
					areasIntersect = false;
			}
		}
	}
	if (areasIntersect) {
		// while testNode is still a inner node (more than one proc), recursively call
		// this method for the children
		if (testNode->_numProcs > 1) {
			getOwningProcs(low, high, decompTree, testNode->_child1, procIDs, neighbHaloAreas);
			getOwningProcs(low, high, decompTree, testNode->_child2, procIDs, neighbHaloAreas);
		}
		else { // leaf node found, add it (and it's area)
			procIDs->push_back(testNode->_owningProc);
			neighbHaloAreas->push_back(testNode->_lowCorner[0] - 1);
			neighbHaloAreas->push_back(testNode->_lowCorner[1] - 1);
			neighbHaloAreas->push_back(testNode->_lowCorner[2] - 1);
			neighbHaloAreas->push_back(testNode->_highCorner[0] + 1);
			neighbHaloAreas->push_back(testNode->_highCorner[1] + 1);
			neighbHaloAreas->push_back(testNode->_highCorner[2] + 1);
		}
	}
}

void KDDecomposition::calcNumParticlesPerCell(ParticleContainer* moleculeContainer) {
	for (int i = 0; i < _globalNumCells * _numParticleTypes; i++)
		_numParticlesPerCell[i] = 0;

	double bBMin[3]; // haloBoundingBoxMin
	for (int dim = 0; dim < 3; dim++) {
		bBMin[dim] = moleculeContainer->getBoundingBoxMin(dim);
	}

	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		for(auto molPtr = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molPtr.isValid(); ++molPtr) {
			int localCellIndex[3]; // 3D Cell index (local)
			int globalCellIdx[3]; // 3D Cell index (global)
			for (int dim = 0; dim < 3; dim++) {
				localCellIndex[dim] = (int) floor((molPtr->r(dim) - bBMin[dim]) / _cellSize[dim]);
				globalCellIdx[dim] = _ownArea->_lowCorner[dim] + localCellIndex[dim];
				if (globalCellIdx[dim] < 0)
					globalCellIdx[dim] += _globalCellsPerDim[dim];
				if (globalCellIdx[dim] >= _globalCellsPerDim[dim])
					globalCellIdx[dim] -= _globalCellsPerDim[dim];
			}
			if(molPtr->componentid() == 0){
				mardyn_assert(static_cast<int>(molPtr->componentid()) <= _numParticleTypes);
				#if defined(_OPENMP)
				#pragma omp atomic
				#endif
				_numParticlesPerCell[_globalCellsPerDim[0] * (globalCellIdx[2] * _globalCellsPerDim[1] + globalCellIdx[1]) + globalCellIdx[0]]++;
			} else {
				#if defined(_OPENMP)
				#pragma omp atomic
				#endif
				_numParticlesPerCell[_globalNumCells + _globalCellsPerDim[0] * (globalCellIdx[2] * _globalCellsPerDim[1] + globalCellIdx[1]) + globalCellIdx[0]]++;
			}
		}
	}
	MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, _numParticlesPerCell.data(), _globalNumCells * _numParticleTypes, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD) );
}

std::vector<int> KDDecomposition::getNeighbourRanks() {
	//global_log->error() << "not implemented \n";
	Simulation::exit(-1);
	return std::vector<int> (0);
}

std::vector<int> KDDecomposition::getNeighbourRanksFullShell() {
	//global_log->error() << "not implemented \n";
	Simulation::exit(-1);
	return std::vector<int> (0);
}

std::vector<CommunicationPartner> KDDecomposition::getNeighboursFromHaloRegion(Domain* domain,
		const HaloRegion& haloRegion, double cutoff) {
//TODO: change this method for support of midpoint rule, half shell, eighth shell, Neutral Territory
	std::vector<CommunicationPartner> communicationPartners;
	int ownLo[DIMgeom];
	int ownHi[DIMgeom];
	double shift[DIMgeom];

	for (unsigned short d = 0; d < DIMgeom; d++) {
		ownLo[d] = _ownArea->_lowCorner[d];
		ownHi[d] = _ownArea->_highCorner[d];
		shift[d] = 0.;
	}

	int regToSendLo[DIMgeom];
	int regToSendHi[DIMgeom];

	double rmin[3], rmax[3];
	for (unsigned int d = 0; d < DIMgeom; d++) {  // put boundary options, such that everything lies within domain
		// TODO: only works for FullShell!
		// (theoretically this could start at offset[d]=-2 or so, if HaloRegion goes over multiple ranks,
		// i.e. if cellsize is less than 1 cutoff radius)
		if (haloRegion.offset[d] < 0 && ownLo[d] == 0) {
			rmin[d] = haloRegion.rmin[d] + domain->getGlobalLength(d);
			rmax[d] = haloRegion.rmax[d] + domain->getGlobalLength(d);
			shift[d] = domain->getGlobalLength(d);
		} else if (haloRegion.offset[d] > 0 && ownHi[d] == _globalCellsPerDim[d] - 1) {
			rmin[d] = haloRegion.rmin[d] - domain->getGlobalLength(d);
			rmax[d] = haloRegion.rmax[d] - domain->getGlobalLength(d);
			shift[d] = -domain->getGlobalLength(d);
		} else {
			rmin[d] = haloRegion.rmin[d];
			rmax[d] = haloRegion.rmax[d];
		}
	}

	getCellIntCoordsFromRegionPeriodic(regToSendLo, regToSendHi, rmin, rmax, domain);

	vector<int> ranks;
	vector<int> ranges;
	_decompTree->getOwningProcs(regToSendLo, regToSendHi, ranks, ranges);
	int numNeighbours = ranks.size();
	auto indexIt = ranges.begin();
	for (int n = 0; n < numNeighbours; ++n) {
		int low[3];
		int high[3];
		bool enlarged[3][2];
		for (int d = 0; d < 3; ++d) {
			low[d] = *(indexIt++);
			high[d] = *(indexIt++);
			if (haloRegion.offset[d] != 0) {
				mardyn_assert(low[d] == high[d]); // TODO: only for FULLSHELL!!!
			}
			enlarged[d][0] = enlarged[d][1] = false;
		}
		for (unsigned int d = 0; d < DIMgeom; d++) { // put boundary options, such that they lie around ownRegion again.
			// TODO: only for FullShell! (theoretically this could start at offset[d]=-2 or so, if HaloRegion goes over multiple ranks
			if (haloRegion.offset[d] < 0 && ownLo[d] == 0) {
				low[d] -= _globalCellsPerDim[d];
				high[d] -= _globalCellsPerDim[d];
			} else if (haloRegion.offset[d] > 0 && ownHi[d] == _globalCellsPerDim[d] - 1) {
				low[d] += _globalCellsPerDim[d];
				high[d] += _globalCellsPerDim[d];
			} else if (haloRegion.offset[d] == 0) {
				if (ownLo[d] < low[d]) {
					low[d]--;
					enlarged[d][0] = true;
				}
				if (ownHi[d] > high[d]) {
					high[d]++;
					enlarged[d][1] = true;
				}
			}
		}

		// region given by the current low-high range is the halo-range
		double haloLow[3];
		double haloHigh[3];
		getCellBorderFromIntCoords(haloLow, haloHigh, low, high);

		for (unsigned int d = 0; d < DIMgeom; d++) {
			// shift for boundaryRegion (one cell)
			low[d] -= haloRegion.offset[d];
			high[d] -= haloRegion.offset[d];
		}

		double boundaryLow[3];
		double boundaryHigh[3];
		getCellBorderFromIntCoords(boundaryLow, boundaryHigh, low, high);

		communicationPartners.emplace_back(ranks[n], haloLow, haloHigh, boundaryLow, boundaryHigh, shift, haloRegion.offset, enlarged);

	}

	return communicationPartners;
}

void KDDecomposition::collectMoleculesInRegion(ParticleContainer* moleculeContainer, const double startRegion[3], const double endRegion[3], vector<Molecule*>& mols) const {
	vector<vector<Molecule*>> threadData;
	vector<int> prefixArray;

	#if defined (_OPENMP)
	#pragma omp parallel shared(mols, threadData)
	#endif
	{
		const int prevNumMols = mols.size();
		const int numThreads = mardyn_get_num_threads();
		const int threadNum = mardyn_get_thread_num();
		auto begin = moleculeContainer->regionIterator(startRegion, endRegion, ParticleIterator::ONLY_INNER_AND_BOUNDARY);

		#if defined (_OPENMP)
		#pragma omp master
		#endif
		{
			threadData.resize(numThreads);
			prefixArray.resize(numThreads + 1);
		}

		#if defined (_OPENMP)
		#pragma omp barrier
		#endif

		for (auto i = begin; i.isValid(); ++i) {
			threadData[threadNum].push_back(new Molecule(*i));
            moleculeContainer->deleteMolecule(i, false); //removeFromContainer = true;
		}

		prefixArray[threadNum + 1] = threadData[threadNum].size();

		#if defined (_OPENMP)
		#pragma omp barrier
		#endif

		//build the prefix array and resize the molecule array
		#if defined (_OPENMP)
		#pragma omp master
		#endif
		{
			int totalNumMols = 0;
			//build the prefix array
			prefixArray[0] = 0;
			for(int i = 1; i <= numThreads; i++){
				prefixArray[i] += prefixArray[i - 1];
				totalNumMols += threadData[i - 1].size();
			}

			//resize the molecule array
			mols.resize(prevNumMols + totalNumMols);
		}

		#if defined (_OPENMP)
		#pragma omp barrier
		#endif

		//reduce the molecules in the molecule array
		int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
		for(int i = 0; i < myThreadMolecules; i++){
			mols[prevNumMols + prefixArray[threadNum] + i] = threadData[threadNum][i];
		}
	}
}

bool KDDecomposition::heteroDecompose(KDNode* fatherNode, KDNode*& ownArea, MPI_Comm commGroup) {
	bool domainTooSmall = false;

	KDNode *bestSubdivision = nullptr;
	domainTooSmall = calculateHeteroSubdivision(fatherNode, bestSubdivision, commGroup);
	double minimalDeviation = FLT_MAX;

	// compute the next subdivision depending on the current rank (either first or second subdivision)
	vector<int> origRanks;
	int newNumProcs;

	if (_rank < bestSubdivision->_child2->_owningProc) {
		origRanks.resize(bestSubdivision->_child1->_numProcs);
		for (int i = 0; i < bestSubdivision->_child1->_numProcs; i++) {
			origRanks[i] = i;  // this group will consist of the first (*iter)->_child1->_numProcs processes/ranks of the current communicator (origGroup)
		}
		newNumProcs = bestSubdivision->_child1->_numProcs;
	} else {
		// ...or the second group (child2) for the MPI communication
		origRanks.resize(bestSubdivision->_child2->_numProcs);
		for (int i = 0; i < bestSubdivision->_child2->_numProcs; i++) {
			origRanks[i] = i + bestSubdivision->_child1->_numProcs; // this group consists of the last (*iter)->_child2->_numProcs processes/ranks of the current communicator (origGroup)
		}
		newNumProcs = bestSubdivision->_child2->_numProcs;
	}

	MPI_Comm newComm;
	MPI_Group origGroup, newGroup;

	MPI_CHECK( MPI_Comm_group(commGroup, &origGroup) );
	MPI_CHECK( MPI_Group_incl(origGroup, newNumProcs, &origRanks[0], &newGroup) );//create new MPI group based on rank (as calculated before)
	MPI_CHECK( MPI_Comm_create(commGroup, newGroup, &newComm) );

	KDNode* newOwnArea = nullptr;
	double deviationChildren[] = {0.0, 0.0};
	if (_rank < bestSubdivision->_child2->_owningProc) {  // compute the subdivision of the first child ...
		// do not use the function call directly in the logical expression, as it may
		// not be executed due to conditional / short-circuit evaluation!
		bool subdomainTooSmall = decompose(bestSubdivision->_child1, newOwnArea, newComm, minimalDeviation);
		deviationChildren[0] = bestSubdivision->_child1->_deviation;
		domainTooSmall = (domainTooSmall || subdomainTooSmall);
	} else {									  // ... or the second child
		mardyn_assert(_rank >= bestSubdivision->_child2->_owningProc);
		bool subdomainTooSmall = decompose(bestSubdivision->_child2, newOwnArea, newComm, minimalDeviation);
		deviationChildren[1] = bestSubdivision->_child2->_deviation;
		domainTooSmall = (domainTooSmall || subdomainTooSmall);
	}
	ownArea = newOwnArea;

	MPI_CHECK( MPI_Group_free(&newGroup));
	MPI_CHECK( MPI_Comm_free(&newComm) );
	MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, deviationChildren, 2, MPI_DOUBLE, _deviationReductionOperation, commGroup));  // reduce the deviations
	bestSubdivision->_child1->_deviation = deviationChildren[0];
	bestSubdivision->_child2->_deviation = deviationChildren[1];
	bestSubdivision->calculateDeviation();


	*fatherNode = *bestSubdivision;  // assignment operator (NOT copy operator) -> also assigns children to fatherNode
	bestSubdivision->_child1 = nullptr;  // remove children from bestSubdivision, otherwise they will be deleted
	bestSubdivision->_child2 = nullptr;  // remove children from bestSubdivision, otherwise they will be deleted
	delete bestSubdivision;

	return domainTooSmall;
}

int KDDecomposition::calculatePartitionRank(){
	int procCount;
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//share the name of all processes between all processes so that everyone can independently determine the partition rank
	char name[MPI_MAX_PROCESSOR_NAME];
	std::fill(name, name+MPI_MAX_PROCESSOR_NAME, ' ');
	int len;
	MPI_Get_processor_name(name, &len);
	std::vector<char> names(MPI_MAX_PROCESSOR_NAME*procCount);
	MPI_Allgather(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
			names.data(), MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD);

	int minDiff = MPI_MAX_PROCESSOR_NAME+1;
	int pos = 0;
	for(int i = 1; i < procCount; ++i){
		char *begin = names.data()+(i-1)*MPI_MAX_PROCESSOR_NAME;
		char *end = names.data()+i*MPI_MAX_PROCESSOR_NAME;
		int diff = std::mismatch(begin, end, end).first - begin;
		std::string test {begin};
		if(diff < minDiff){
			minDiff = diff;
			pos = i;
		}
	}
	return pos;
}

bool KDDecomposition::calculateHeteroSubdivision(KDNode* node, KDNode*& optimalNode, MPI_Comm commGroup) {
	bool domainTooSmall = false;
	vector<vector<double> > costsLeft(3);
	vector<vector<double> > costsRight(3);
	{
		MPI_Group group1, group2; //one group for each partition
		std::vector<int> ranks1 {};
		for(int i = 0; i < _partitionRank; ++i){
			ranks1.push_back(i);
		}
		mardyn_assert(static_cast<int>(ranks1.size()) == _partitionRank);
		MPI_Comm newComm1, newComm2;
		MPI_Group newGroup1, newGroup2;
		MPI_Group origGroup;
		int rank;
		MPI_CHECK( MPI_Comm_rank(commGroup, &rank));
		MPI_CHECK( MPI_Comm_group(commGroup, &origGroup) );
		MPI_CHECK( MPI_Group_incl(origGroup, ranks1.size(), ranks1.data(), &newGroup1) );//create new MPI group based on rank (as calculated before)
		MPI_CHECK( MPI_Comm_create(commGroup, newGroup1, &newComm1) );

		MPI_CHECK( MPI_Group_excl(origGroup, ranks1.size(), ranks1.data(), &newGroup2));//create new MPI group based on rank (as calculated before)
		MPI_CHECK( MPI_Comm_create(commGroup, newGroup2, &newComm2) );

		int currNumProcs = node->_numProcs;
		vector<vector<double>> dummyCosts {3};
		//calculate the load for one of the sides of each splitting plane
		if(rank < _partitionRank){
			node->_numProcs = _partitionRank;
			calculateCostsPar(node, costsLeft, dummyCosts, newComm1);
		} else {
			node->_numProcs = _numProcs - _partitionRank;
			calculateCostsPar(node, dummyCosts, costsRight, newComm2);
		}
		node->_numProcs = currNumProcs;

		//just bring not calculated costs to the right size
		if(rank < _partitionRank){
			costsRight = std::move(dummyCosts);
		} else {
			costsLeft = std::move(dummyCosts);
		}

		/*
		 * Communicate the measured values between all involved processes
		 */
		for(std::vector<double>& vec : costsLeft){
			MPI_Bcast(vec.data(), vec.size(), MPI_DOUBLE, 0, commGroup);

		}

		for(std::vector<double>& vec : costsRight){
			MPI_Bcast(vec.data(), vec.size(), MPI_DOUBLE, _partitionRank, commGroup);
		}

		node->_numProcs = currNumProcs;
	}

	//the ratio is the number of processes in one partition divided by the number of procs in the other one
	double leftRightLoadRatio =  double(_partitionRank)/(node->_numProcs - _partitionRank);

	size_t max = costsLeft[0].size();
	size_t maxInd = 0;
	for (unsigned int dim = 1; dim < 3; dim++){
		if (costsLeft[dim].size() > max){
			max = costsLeft[dim].size();
			maxInd = dim;
		}
	}
	size_t biggestDim = maxInd;

	if (costsLeft[biggestDim].size()<=2){
		global_log->error_always_output() << "The domain is far to small!";
		Simulation::exit(1);
	}

	int startIndex = 1;
	int maxEndIndex = node->_highCorner[biggestDim] - node->_lowCorner[biggestDim] - 1;
	int endIndex = maxEndIndex;

	double minError = fabs(costsLeft[biggestDim][0] / costsRight[biggestDim][0] - leftRightLoadRatio);
	size_t index = 0;
	double error;
	for (size_t i = 1; i < costsLeft[biggestDim].size(); ++i) {
		auto left = costsLeft[biggestDim][i];
		auto right = costsRight[biggestDim][i];

		error = fabs(costsLeft[biggestDim][i] / costsRight[biggestDim][i] - leftRightLoadRatio);
		if (error < minError) {
			minError = error;
			index = i;
		}
	}

	startIndex = std::max((size_t)startIndex, index);
	endIndex = std::min(startIndex + 1, endIndex);
	startIndex = std::min(endIndex-1, startIndex);

	global_log->debug() << "splitLoad: startindex " << index << " of " << costsLeft[biggestDim].size() << std::endl;


	const int i = startIndex;

	double optCostPerProc = (costsLeft[biggestDim][i] + costsRight[biggestDim][i]) / ((double) node->_numProcs);
	int optNumProcsLeft;

	int numProcsLeft = _partitionRank;

	optimalNode = new KDNode(*node);
	if (optimalNode->_level == 0) {
		optimalNode->_optimalLoadPerProcess = optCostPerProc;
	}

	optimalNode->split(biggestDim, node->_lowCorner[biggestDim] + i, numProcsLeft);
	if ( (unsigned int) (optimalNode->_child1->_numProcs + optimalNode->_child2->_numProcs) >
	        (optimalNode->_child1->getNumMaxProcs() + optimalNode->_child2->getNumMaxProcs())) {
		global_log->error() << "Domain is not resolvable at all!" << std::endl;
		Simulation::exit(1);
	}

	while ( (! optimalNode->_child1->isResolvable()) && optimalNode->_child2->isResolvable()) {
		// shift procs to child 2, adapt owner of child 2
		optimalNode->_child1->_numProcs--;
		optimalNode->_child2->_numProcs++;
		optimalNode->_child2->_owningProc--;
		optimalNode->_child2->_nodeID = optimalNode->_nodeID + 2 * optimalNode->_child1->_numProcs;
		domainTooSmall = true;
	}

	while ( optimalNode->_child1->isResolvable() && (! optimalNode->_child2->isResolvable())) {
		// shift procs to child 1, , adapt owner of child 2
		optimalNode->_child1->_numProcs++;
		optimalNode->_child2->_numProcs--;
		optimalNode->_child2->_owningProc++;
		optimalNode->_child2->_nodeID = optimalNode->_nodeID + 2 * optimalNode->_child1->_numProcs;
		domainTooSmall = true;
	}

	// In this place, MBu had some processor shifting in his algorithm.
	// I believe it to be unnecessary, as the ratio of left and right processors
	// is chosen according to the load ratio (I use round instead of floor).
	if ((optimalNode->_child1->_numProcs <= 0 || optimalNode->_child1->_numProcs >= node->_numProcs) ||
			(optimalNode->_child2->_numProcs <= 0 || optimalNode->_child2->_numProcs >= node->_numProcs) ){
		//continue;
		global_log->error_always_output() << "ERROR in calculateHeteroSubdivision(), part of the domain was not assigned to a proc" << endl;
		Simulation::exit(1);
	}
	mardyn_assert( optimalNode->_child1->isResolvable() && optimalNode->_child2->isResolvable() );

	optimalNode->_child1->_load = costsLeft[biggestDim][i];
	optimalNode->_child2->_load = costsRight[biggestDim][i];
	optimalNode->_load = costsLeft[biggestDim][i] + costsRight[biggestDim][i];
	optimalNode->calculateDeviationLowerBound();

	return domainTooSmall;
}
void KDDecomposition::printTree(std::ostream& ostream) {
	_decompTree->printTree("", ostream);
}
