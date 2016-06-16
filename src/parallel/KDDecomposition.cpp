#include "KDDecomposition.h"

#include <cfloat>
#include <sstream>
#include <fstream>
#include <climits>

#include "Domain.h"
#include "KDNode.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "particleContainer/ParticleContainer.h"
#include "ParticleData.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"
#include "particleContainer/adapter/FlopCounter.h"

#include <cmath>

using namespace std;
using Log::global_log;

//#define DEBUG_DECOMP

KDDecomposition::KDDecomposition() :
		_globalNumCells(1), _decompTree(NULL), _ownArea(NULL), _numParticlesPerCell(NULL), _steps(0), _frequency(1.), _cutoffRadius(
				1.), _fullSearchThreshold(8), _totalMeanProcessorSpeed(1.), _totalProcessorSpeed(1.), _processorSpeedUpdateCount(0),
				_heterogeneousSystems(false), _splitBiggest(true), _forceRatio(false){
	bool before = global_log->get_do_output();
	global_log->set_mpi_output_all();
	global_log->debug() << "KDDecomposition: Rank " << _rank << " executing file " << global_simulation->getName() << std::endl;
	global_log->set_do_output(before);
}

KDDecomposition::KDDecomposition(double cutoffRadius, Domain* domain, int updateFrequency, int fullSearchThreshold, bool hetero, bool cutsmaller, bool forceRatio) :
		_steps(0), _frequency(updateFrequency), _fullSearchThreshold(fullSearchThreshold), _totalMeanProcessorSpeed(1.),
		_totalProcessorSpeed(1.), _processorSpeedUpdateCount(0), _heterogeneousSystems(hetero), _splitBiggest(!cutsmaller), _forceRatio(forceRatio) {
	bool before = global_log->get_do_output();
	global_log->set_mpi_output_all();
	global_log->debug() << "KDDecomposition: Rank " << _rank << " executing file " << global_simulation->getName() << std::endl;
	global_log->set_do_output(before);
	_cutoffRadius = cutoffRadius;

	int lowCorner[KDDIM] = {0};
	int highCorner[KDDIM] = {0};
	bool coversWholeDomain[KDDIM];
	_globalNumCells = 1;
	
	for (int dim = 0; dim < KDDIM; dim++) {
		_globalCellsPerDim[dim] = (int) floor(domain->getGlobalLength(dim) / cutoffRadius);
		_globalNumCells *= _globalCellsPerDim[dim];
		highCorner[dim] = _globalCellsPerDim[dim] - 1;
		_cellSize[dim] = domain->getGlobalLength(dim) / ((double) _globalCellsPerDim[dim]);
		coversWholeDomain[dim] = true;
	}
	
	_numParticlesPerCell = new unsigned int[_globalNumCells];
	for (int i = 0; i < _globalNumCells; i++)
		_numParticlesPerCell[i] = 0;

	// create initial decomposition
	// ensure that enough cells for the number of procs are available
	_decompTree = new KDNode(_numProcs, lowCorner, highCorner, 0, 0, coversWholeDomain, 0);
	if (!_decompTree->isResolvable()) {
		global_log->error() << "KDDecompsition not possible. Each process needs at least 8 cells." << endl;
		global_log->error() << "The number of Cells is only sufficient for " << _decompTree->getNumMaxProcs() << " Procs!" << endl;
		barrier(); // the messages above are only promoted to std::out if we have the barrier somehow...
		global_simulation->exit(-1);
	}
	_decompTree->buildKDTree();
	_ownArea = _decompTree->findAreaForProcess(_rank);

	// initialize the mpi data type for particles once in the beginning
	KDNode::initMPIDataType();

	global_log->info() << "Created KDDecomposition with updateFrequency=" << _frequency << ", fullSearchThreshold=" << _fullSearchThreshold << endl;

#ifdef DEBUG_DECOMP
	global_log->info() << "Initial Decomposition: " << endl;
	if (_rank == 0) {
		_decompTree->printTree("");
	}
#endif
}

KDDecomposition::~KDDecomposition() {
	delete[] _numParticlesPerCell;
//	_decompTree->serialize(string("kddecomp.dat"));
	if (_rank == 0) {
		_decompTree->plotNode("kddecomp.vtu", &_processorSpeeds);
	}
	delete _decompTree;
	KDNode::shutdownMPIDataType();
}

void KDDecomposition::readXML(XMLfileUnits& xmlconfig) {
	/* TODO: Maybe add decomposition dimensions, default auto. */
	xmlconfig.getNodeValue("updateFrequency", _frequency);
	global_log->info() << "KDDecomposition update frequency: " << _frequency << endl;
	xmlconfig.getNodeValue("fullSearchThreshold", _fullSearchThreshold);
	global_log->info() << "KDDecomposition full search threshold: " << _fullSearchThreshold << endl;
	xmlconfig.getNodeValue("heterogeneousSystems", _heterogeneousSystems);
	global_log->info() << "KDDecomposition for heterogeneous systems?: " << (_heterogeneousSystems?"yes":"no") << endl;
	xmlconfig.getNodeValue("splitBiggestDimension", _splitBiggest);
	global_log->info() << "KDDecomposition splits along biggest domain?: " << (_splitBiggest?"yes":"no") << endl;
	xmlconfig.getNodeValue("forceRatio", _forceRatio);
	global_log->info() << "KDDecomposition forces load/performance ratio?: " << (_forceRatio?"yes":"no") << endl;
}

int KDDecomposition::getNonBlockingStageCount(){
	return 3;
}

void KDDecomposition::prepareNonBlockingStage(bool /*forceRebalancing*/,
		ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber) {
	const bool removeRecvDuplicates = true;
	DomainDecompMPIBase::prepareNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_AND_HALO_COPIES, removeRecvDuplicates);
}

void KDDecomposition::finishNonBlockingStage(bool /*forceRebalancing*/,
		ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber) {
	const bool removeRecvDuplicates = true;
	DomainDecompMPIBase::finishNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_AND_HALO_COPIES, removeRecvDuplicates);
}

//check whether or not to do rebalancing in the specified step
bool doRebalancing(bool forceRebalancing, size_t steps, int frequency){
	return forceRebalancing or steps % frequency == 0 or steps <= 1;
}

bool KDDecomposition::queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* /*moleculeContainer*/, Domain* /*domain*/){
	return not doRebalancing(forceRebalancing, _steps, _frequency);
}

void KDDecomposition::balanceAndExchange(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain) {
	const bool rebalance = doRebalancing(forceRebalancing, _steps, _frequency);
	_steps++;
	const bool removeRecvDuplicates = true;

	if (rebalance == false) {
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES, removeRecvDuplicates);
	} else {
		global_log->info() << "KDDecomposition: rebalancing..." << endl;

		if (_steps != 1) {
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY, removeRecvDuplicates);
		}
		moleculeContainer->deleteOuterParticles();

		KDNode * newDecompRoot = NULL;
		KDNode * newOwnLeaf = NULL;

		getNumParticles(moleculeContainer);
		constructNewTree(newDecompRoot, newOwnLeaf, moleculeContainer);
		bool migrationSuccessful = migrateParticles(*newDecompRoot, *newOwnLeaf, moleculeContainer);
		if (not migrationSuccessful) {
			global_log->error() << "A problem occurred during particle migration between old decomposition and new decomposition of the KDDecomposition." << endl;
			global_log->error() << "Aborting. Please save your input files and last available checkpoint and contact TUM SCCS." << endl;
			global_simulation->exit(1);
		}
		delete _decompTree;
		_decompTree = newDecompRoot;
//		delete _ownArea; dont delete! this is a pointer only to one of the objects in the whole tree, not a real object
		_ownArea = newOwnLeaf;
		initCommunicationPartners(_cutoffRadius, domain);

		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES, removeRecvDuplicates);
	}
}

void KDDecomposition::initCommunicationPartners(double cutoffRadius, Domain * domain) {

//	if(_neighboursInitialized) {
//		return;
//	}
//	_neighboursInitialized = true;

	int ownLo[DIM];
	int ownHi[DIM];

	for (unsigned short d = 0; d < DIM; d++) {
		ownLo[d] = _ownArea->_lowCorner[d];
		ownHi[d] = _ownArea->_highCorner[d];

		_neighbours[d].clear();
	}

	for (unsigned short dimension = 0; dimension < DIM; dimension++) {
		if(_coversWholeDomain[dimension]) {
			// nothing to do;
			continue;
		}

		for (int direction = LOWER; direction <= HIGHER; direction++) {
			double shift[DIM];
			for (int i = 0; i < 3; ++i) {
				shift[i] = 0.0;
			}

			int regToSendLo[DIM];
			int regToSendHi[DIM];

			for (int i = 0; i < DIM; ++i) {
				regToSendLo[i] = ownLo[i];
				regToSendHi[i] = ownHi[i];
			}

			if (ownLo[dimension] == 0) {
				regToSendLo[dimension] = _globalCellsPerDim[dimension];
			} else if (ownHi[dimension] == _globalCellsPerDim[dimension] - 1) {
				regToSendHi[dimension] = -1;
			}

			switch (direction) {
			case LOWER:
				--regToSendLo[dimension];
				regToSendHi[dimension] = regToSendLo[dimension];
				if (ownLo[dimension] == 0) {
					shift[dimension] = domain->getGlobalLength(dimension);
				}
				break;
			case HIGHER:
				++regToSendHi[dimension];
				regToSendLo[dimension] = regToSendHi[dimension];
				if (ownHi[dimension] == _globalCellsPerDim[dimension] - 1) {
					shift[dimension] = -domain->getGlobalLength(dimension);
				}
				break;
			}

			vector<int> ranks;
			vector<int> ranges;
			_decompTree->getOwningProcs(regToSendLo, regToSendHi, ranks, ranges);
			int numNeighbours = ranks.size();
			vector<int>::iterator indexIt = ranges.begin();
			for (int n = 0; n < numNeighbours; ++n) {
				int low[3];
				int high[3];
				for (int d=0; d < 3; ++d) {
					low[d] = *(indexIt++);
					high[d] = *(indexIt++);
					if (d == dimension) {
						assert(low[d] == high[d]);
					}
				}
				switch (direction) {
				case LOWER:
					if (low[dimension] + 1 != ownLo[dimension]) {
						low[dimension] = high[dimension] = ownLo[dimension] - 1;
					}
					break;
				case HIGHER:
					if (high[dimension] - 1 != ownHi[dimension]) {
						high[dimension] = low[dimension] = ownHi[dimension] + 1;
					}
					break;
				}

				// enlarge in two "other" dimensions
				for (unsigned short d=0; d < 3; ++d) {
					if(d != dimension) {
						--low[d];
						++high[d];
					}
				}

				// region given by the current low-high range is the halo-range
				double haloLow[3];
				double haloHigh[3];
				getCellBorderFromIntCoords(haloLow, haloHigh, low, high);

				switch (direction) {
				case LOWER:
					low[dimension]++;
					high[dimension]++;
					break;
				case HIGHER:
					low[dimension]--;
					high[dimension]--;
					break;
				}

				double boundaryLow[3];
				double boundaryHigh[3];
				getCellBorderFromIntCoords(boundaryLow, boundaryHigh, low, high);

				_neighbours[dimension].push_back(CommunicationPartner(ranks[n], haloLow, haloHigh, boundaryLow, boundaryHigh, shift));

			}

		}
	}
}

void KDDecomposition::getCellBorderFromIntCoords(double * lC, double * hC, int lo[3], int hi[3]) const {
	for(int d = 0 ; d < 3; ++d) {
		lC[d] = lo[d] * _cellSize[d];
		hC[d] = (hi[d] + 1) * _cellSize[d];
	}
}

bool KDDecomposition::migrateParticles(const KDNode& newRoot, const KDNode& newOwnLeaf, ParticleContainer* moleculeContainer) const {
	// 1. compute which processes we will receive from
	// 2. issue Irecv calls
	// 3. compute which prcesses we will send to
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
		vector<int>::iterator indexIt = indices.begin();
		numProcsRecv = ranks.size(); // value may change from ranks.size(), see "numProcsSend--" below
		recvPartners.reserve(numProcsRecv);
		for (unsigned i = 0; i < ranks.size(); ++i) {
			int partnerRank = ranks[i];

			if (partnerRank != _rank) {
				recvPartners.push_back(CommunicationPartner(partnerRank));
			}

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
					}
				}
			}

			if (partnerRank != _rank) {
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

		vector<int>::iterator indexIt = indices.begin();
		numProcsSend = ranks.size(); // value may change from ranks.size(), see "numProcsSend--" below
		sendPartners.reserve(numProcsSend);
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
				sendPartners.push_back(CommunicationPartner(partnerRank, leavingLow, leavingHigh));
				const bool removeFromContainer = true;
				sendPartners.back().initSend(moleculeContainer, _comm, _mpiParticleType, LEAVING_ONLY, removeFromContainer); // molecules have been taken out of container
			} else {
				moleculeContainer->getRegionSimple(leavingLow, leavingHigh, migrateToSelf, true);
				// decrement numProcsSend for further uses:
				assert(willMigrateToSelf == true);
				numProcsSend--;
			}
		}
	}
	assert(moleculeContainer->getNumberOfParticles() == 0ul);
	double newBoxMin[3];
	double newBoxMax[3];
	for (int dim = 0; dim < 3; dim++) {
		newBoxMin[dim] = (newOwnLeaf._lowCorner[dim]) * _cellSize[dim];
		newBoxMax[dim] = (newOwnLeaf._highCorner[dim] + 1) * _cellSize[dim];
	}
	moleculeContainer->rebuild(newBoxMin, newBoxMax);

	global_log->set_mpi_output_all();
	double waitCounter = 1.0;
	double deadlockTimeOut = 60.0;
	bool allDone = false;
	double startTime = MPI_Wtime();
	bool migrateToSelfDone = not willMigrateToSelf;

	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (int i = 0; i < numProcsSend; ++i) {
			allDone &= sendPartners[i].testSend();
		}

		if (migrateToSelfDone != true) {
			const int numMolsMigToSelf = migrateToSelf.size();
			for (int i = 0; i < numMolsMigToSelf; i++) {
				moleculeContainer->addParticlePointer(migrateToSelf[i], false, false);
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
			global_log->warning() << "Deadlock warning: Rank " << _rank
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
			global_log->warning() << "Deadlock error: Rank " << _rank
					<< " is waiting for more than " << deadlockTimeOut
					<< " seconds" << std::endl;
			for (int i = 0; i < numProcsSend; ++i) {
				sendPartners[i].deadlockDiagnosticSend();
			}
			for (int i = 0; i < numProcsRecv; ++i) {
				recvPartners[i].deadlockDiagnosticRecv();
			}
			global_log->warning() << "aborting" << std::endl;
			break;
		}

	} // while not allDone

	moleculeContainer->update();

	global_log->set_mpi_output_root(0);

	int isOK = allDone;

	MPI_Allreduce(MPI_IN_PLACE, &isOK, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	bool success = false;

	if(isOK == _numProcs) {
		success = true;
	}

	return success;
}

void KDDecomposition::constructNewTree(KDNode *& newRoot, KDNode *& newOwnLeaf, ParticleContainer* moleculeContainer) {

	newRoot = new KDNode(_numProcs, &(_decompTree->_lowCorner[0]), &(_decompTree->_highCorner[0]), 0, 0, _decompTree->_coversWholeDomain, 0);
	KDNode * toCleanUp = newRoot;

	updateMeanProcessorSpeeds(_processorSpeeds,_accumulatedProcessorSpeeds, moleculeContainer);

	if (decompose(newRoot, newOwnLeaf, MPI_COMM_WORLD)) {
		global_log->warning() << "Domain too small to achieve a perfect load balancing" << endl;
	}

	completeTreeInfo(newRoot, newOwnLeaf);
	delete toCleanUp;
	for (int d = 0; d < 3; ++d) {
		_coversWholeDomain[d] = newOwnLeaf->_coversWholeDomain[d];
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
		newRoot->printTree("");
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
		if (processorSpeeds.size() == 0) {
			processorSpeeds.resize(_numProcs);
			accumulatedProcessorSpeeds.clear();
		}
		for (int i = 0; i < _numProcs; i++) {
			processorSpeeds[i] = collCommGetDouble();
		}
		collCommFinalize();
	}
	else{
		if (processorSpeeds.size() == 0) {
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
		assert(processorSpeeds[i] > 0);
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

void KDDecomposition::printDecomp(string filename, Domain* domain) {
	if (_rank == 0) {
		ofstream povcfgstrm(filename.c_str());
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " " << domain->getGlobalLength(2) << endl;
		povcfgstrm << "decompData Regions" << endl;
		povcfgstrm.close();
	}

	for (int process = 0; process < _numProcs; process++) {
		if (_rank == process) {
			ofstream povcfgstrm(filename.c_str(), ios::app);
			povcfgstrm << getBoundingBoxMin(0,domain) << " " << getBoundingBoxMin(1,domain) << " "
			           << getBoundingBoxMin(2,domain) << " " << getBoundingBoxMax(0,domain) << " "
			           << getBoundingBoxMax(1,domain) << " " << getBoundingBoxMax(2,domain) << endl;
			povcfgstrm.close();
		}
		barrier();
	}
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
		ptrToAllNodes[nodeID]->_expectedDeviation = mpiKDNode.getExpectedDeviation();
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
	filestream << " * " << "load=" << node->_load << " optLoad=" << node->_optimalLoadPerProcess << " expDev=" << node->_expectedDeviation << " minDev=" << minDev << endl;
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
		assert(fatherNode->_owningProc == _rank);
		ownArea = fatherNode;
		fatherNode->calculateDeviation(&_processorSpeeds, _totalMeanProcessorSpeed);
		return domainTooSmall;
	}

	std::list<KDNode*> subdivisions;
	domainTooSmall = calculateAllSubdivisions(fatherNode, subdivisions, commGroup);
	assert(subdivisions.size() > 0);

	KDNode* bestSubdivision = NULL;
	double minimalDeviation = globalMinimalDeviation;
	list<KDNode*>::iterator iter = subdivisions.begin();
	int iterations = 0;
	int log2proc = 0;// calculates the logarithm of _numProcs (base 2)
	while ((_numProcs >> log2proc) > 1) {
		log2proc++;
	}
	// if we are near the root of the tree, we just take the first best subdivision
	// this is even valid for heterogeneous balancing, since the subdivisions are calculated accordingly.
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
               "level=" << fatherNode->_level << " #divisions=" << subdivisions.size() << endl;
#endif

	while (iter !=  subdivisions.end() && (iterations < maxIterations) && (*iter)->_expectedDeviation < minimalDeviation) {
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

		KDNode* newOwnArea = NULL;
		double deviationChildren[] = {0.0, 0.0};

		if (_rank < (*iter)->_child2->_owningProc) {  // compute the subdivision of the first child ...
			// do not use the function call directly in the logical expression, as it may
			// not be executed due to conditional / short-circuit evaluation!
			bool subdomainTooSmall = decompose((*iter)->_child1, newOwnArea, newComm, minimalDeviation);
			deviationChildren[0] = (*iter)->_child1->_deviation;
			domainTooSmall = (domainTooSmall || subdomainTooSmall);
		} else {									  // ... or the second child
			assert(_rank >= (*iter)->_child2->_owningProc);
			bool subdomainTooSmall = decompose((*iter)->_child2, newOwnArea, newComm, minimalDeviation);
			deviationChildren[1] = (*iter)->_child2->_deviation;
			domainTooSmall = (domainTooSmall || subdomainTooSmall);
		}

		MPI_CHECK( MPI_Group_free(&newGroup));
		MPI_CHECK( MPI_Comm_free(&newComm) );
		MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, deviationChildren, 2, MPI_DOUBLE, MPI_SUM, commGroup));  // reduce the deviations
		// TODO hetero: check whether there really has to be an allreduce (sum), since each of the processes should already possess the deviations of all its children.
		//       with the allreduce (sum) implementation trees with a low maximal level are preferred (balanced trees).
		//		 shouldn't this better be an average????
		(*iter)->_child1->_deviation = deviationChildren[0];
		(*iter)->_child2->_deviation = deviationChildren[1];
		(*iter)->calculateDeviation();

#ifdef DEBUG_DECOMP
		for (int i = 0; i < fatherNode->_level; i++) { filestream << "   ";}
		filestream << "   deviation=" << (*iter)->_deviation << " (ch1:" << deviationChildren[0] << "ch2:" << deviationChildren[1] << endl;
#endif
		if ((*iter)->_deviation < minimalDeviation) {
			delete bestSubdivision;// (deleting of NULL is ok and does not produce errors, since delete checks for NULL)
			bestSubdivision = *iter;
			minimalDeviation = (*iter)->_deviation;
			ownArea = newOwnArea;
		} else {
			delete *iter;
		}
		iter++;
	}

	while (iter !=  subdivisions.end()) {
		delete *iter;
		iter++;
	}

	// reassign children and delete cloned node, if a better solution
	// was found in this subtree.
	if (bestSubdivision == NULL) {
		fatherNode->_deviation = FLT_MAX;
	} else {
		*fatherNode = *bestSubdivision;
		bestSubdivision->_child1 = NULL;
		bestSubdivision->_child2 = NULL;
		delete bestSubdivision;
	}
#ifdef DEBUG_DECOMP
	filestream.close();
#endif
	return domainTooSmall;
}


bool KDDecomposition::calculateAllSubdivisions(KDNode* node, std::list<KDNode*>& subdividedNodes, MPI_Comm commGroup) {
	bool domainTooSmall = false;
	vector<vector<double> > costsLeft(3);
	vector<vector<double> > costsRight(3);
	calculateCostsPar(node, costsLeft, costsRight, commGroup);
	global_log->debug() << "calculateAllSubdivisions: " << std::endl;
	double leftRightLoadRatio = 1.;  // divide load 50/50 -> this can be changed later on if the topology of the system should be taken into account.
	// update leftRightRatio to something meaningful (that is representable by the compute power distribution of the processes:

	// speed that is close to LeftRightRatio
	double optimalSpeed = (_accumulatedProcessorSpeeds[node->_owningProc + node->_numProcs]
			- _accumulatedProcessorSpeeds[node->_owningProc]) * leftRightLoadRatio / (1. + leftRightLoadRatio);
	double searchSpeed = optimalSpeed + _accumulatedProcessorSpeeds[node->_owningProc];
	std::vector<double>::iterator iter = std::lower_bound(_accumulatedProcessorSpeeds.begin() + node->_owningProc,
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



	bool splitLoad = true;  // indicates, whether to split the domain according to the load
							// or whether the domain should simply be split in half and the number of processes should be distributed accordingly.

	size_t dimInit = 0;
	size_t dimEnd = 3;
	if(_splitBiggest){
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
		if (costsLeft[dim].size()==0){
			continue;
		}
		// if a node has some more processors, we probably don't have to find the
		// "best" partitioning, but it's sufficient to divide the node in the middle
		// and shift the processor count according to the load imbalance.

		// loop only from 1 to max-1 (instead 0 to max) to avoid 1-cell-regions
		int startIndex = 1;
		int maxEndIndex = node->_highCorner[dim] - node->_lowCorner[dim] - 1;
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
				startIndex = min(endIndex-1, startIndex);
				global_log->debug() << "splitLoad: startindex " << index << " of " << costsLeft[dim].size() <<std::endl;
			} else {  // If we have more than _fullSearchThreshold processes left, we split the domain in half.
				startIndex = max(startIndex, (node->_highCorner[dim] - node->_lowCorner[dim] - 1) / 2);
				endIndex = min(endIndex, startIndex + 1);
			}
		}


		int i = startIndex;
		while ( ((i < endIndex) || (subdividedNodes.size() == 0 && i < maxEndIndex))) {
			//for (int i = startIndex; i < endIndex; i++) {

			double optCostPerProc = (costsLeft[dim][i] + costsRight[dim][i]) / ((double) node->_numProcs);
			int optNumProcsLeft;
			if (splitLoad) {  // if we split the load in a specific ratio, numProcsLeft is calculated differently
				if(_accumulatedProcessorSpeeds.size()==0){
					global_log->error() << "no processor speeds given" << std::endl;
					global_simulation->exit(-1);
				}
				double optimalLoad = (_accumulatedProcessorSpeeds[node->_owningProc + node->_numProcs] - _accumulatedProcessorSpeeds[node->_owningProc]) * leftRightLoadRatio
						/ (1. + leftRightLoadRatio);
				double searchLoad = optimalLoad + _accumulatedProcessorSpeeds[node->_owningProc];
				std::vector<double>::iterator iter = std::lower_bound(_accumulatedProcessorSpeeds.begin() + node->_owningProc,
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

			KDNode* clone = new KDNode(*node);
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

			// In this place, MBu had some processor shifting in his algorithm.
			// I believe it to be unnecessary, as the ratio of left and right processors
			// is chosen according to the load ratio (I use round instead of floor).
			if ((clone->_child1->_numProcs <= 0 || clone->_child1->_numProcs >= node->_numProcs) ||
					(clone->_child2->_numProcs <= 0 || clone->_child2->_numProcs >= node->_numProcs) ){
				//continue;
				global_log->error_always_output() << "ERROR in calculateAllSubdivisions(), part of the domain was not assigned to a proc" << endl;
				global_simulation->exit(1);
			}
			assert( clone->_child1->isResolvable() && clone->_child2->isResolvable() );

			clone->_child1->_load = costsLeft[dim][i];
			clone->_child2->_load = costsRight[dim][i];
			clone->_load = costsLeft[dim][i] + costsRight[dim][i];
			clone->calculateExpectedDeviation(&_accumulatedProcessorSpeeds);

			// sort node according to expected deviation
			list<KDNode*>::iterator iter = subdividedNodes.begin();
			while (iter != subdividedNodes.end() && ((*iter)->_expectedDeviation < clone->_expectedDeviation)) {
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
					int numParts = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, i_dim, i_dim1, i_dim2, area)];

					// #######################
					// ## Cell Costs        ##
					// #######################
					cellCosts[dim][i_dim] += (double) (numParts * numParts);

					// all Neighbours
					for (int neighInd_divDim = i_dim - 1; neighInd_divDim <= i_dim + 1; neighInd_divDim++) {
						int nI_dim = neighInd_divDim;
						if (nI_dim + area->_lowCorner[dim] < 0)
							nI_dim = _globalCellsPerDim[dim] - 1;
						if (nI_dim + area->_lowCorner[dim] >= _globalCellsPerDim[dim])
							nI_dim = 0;

						for (int neighInd_dim1 = i_dim1 - 1; neighInd_dim1 <= i_dim1 + 1; neighInd_dim1++) {
							// adjust index in case of periodic boundary
							int nI_dim1 = neighInd_dim1;
							if (nI_dim1 + area->_lowCorner[dim1] < 0)
								nI_dim1 = _globalCellsPerDim[dim1] - 1;
							if (nI_dim1 + area->_lowCorner[dim1] >= _globalCellsPerDim[dim1])
								nI_dim1 = 0;

							for (int neighInd_dim2 = i_dim2 - 1; neighInd_dim2 <= i_dim2 + 1; neighInd_dim2++) {
								// adjust index in case of periodic boundary
								int nI_dim2 = neighInd_dim2;
								if (nI_dim2 + area->_lowCorner[dim2] < 0)
									nI_dim2 = _globalCellsPerDim[dim2] - 1;
								if (nI_dim2 + area->_lowCorner[dim2] >= _globalCellsPerDim[dim2])
									nI_dim2 = 0;
								// count only forward neighbours
								int numPartsNeigh = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, nI_dim, nI_dim1, nI_dim2, area)];
								cellCosts[dim][i_dim] += 0.5 * (double) (numParts * numPartsNeigh);
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
		if (not (coversWholeDomain[dim])) { // TODO: check wheter this can be removed || _ownArea->_coversWholeDomain[dim])){
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

void KDDecomposition::getNumParticles(ParticleContainer* moleculeContainer) {
	for (int i = 0; i < _globalNumCells; i++)
		_numParticlesPerCell[i] = 0;

	int count = 0;
	double bBMin[3]; // haloBoundingBoxMin
    /* TODO: We do not use values form bBMax anywhere ... */
	//double bBMax[3]; // haloBoundingBoxMax
	for (int dim = 0; dim < 3; dim++) {
		bBMin[dim] = moleculeContainer->getBoundingBoxMin(dim);// - moleculeContainer->get_halo_L(dim);
		//bBMax[dim] = moleculeContainer->getBoundingBoxMax(dim);// + moleculeContainer->get_halo_L(dim);
	}
	Molecule* molPtr = moleculeContainer->begin();
	while (molPtr != moleculeContainer->end()) {
		int cellIndex[3]; // 3D Cell index (local)
		int globalCellIdx[3]; // 3D Cell index (global)

		for (int dim = 0; dim < 3; dim++) {
			cellIndex[dim] = (int) floor((molPtr->r(dim) - bBMin[dim]) / _cellSize[dim]);
			globalCellIdx[dim] = _ownArea->_lowCorner[dim] + cellIndex[dim];
			if (globalCellIdx[dim] < 0)
				globalCellIdx[dim] += _globalCellsPerDim[dim];
			if (globalCellIdx[dim] >= _globalCellsPerDim[dim])
				globalCellIdx[dim] -= _globalCellsPerDim[dim];
		}

		_numParticlesPerCell[_globalCellsPerDim[0] * (globalCellIdx[2] * _globalCellsPerDim[1] + globalCellIdx[1]) + globalCellIdx[0]]++;
		molPtr = moleculeContainer->next();
		count++;
	}
	MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, _numParticlesPerCell, _globalNumCells, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD) );

}

