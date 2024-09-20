/*
 * UniformPseudoParticleContainer.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: tchipevn
 */

#include "UniformPseudoParticleContainer.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"
#include "bhfmm/HaloBufferNoOverlap.h"
#include "bhfmm/HaloBufferOverlap.h"
#include <string>
#include <algorithm>
#include <array>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <bhfmm/FastMultipoleMethod.h>

namespace bhfmm {

#ifndef WIGNER
#define WIGNER 0 // 0: original, 1: Wigner rotation acceleration
#endif

#define HALOSIZE 2(_globalLevel == 1 and _fuseGlobalCommunication)

#define IsOdd(x) ((x) & 1)
#define ToEven(x) ((x) & ~1)

UniformPseudoParticleContainer::UniformPseudoParticleContainer(
		double domainLength[3],
		double bBoxMin[3],
		double bBoxMax[3],
		double LJCellLength[3],
		unsigned LJSubdivisionFactor,
		int orderOfExpansions,
		ParticleContainer *ljContainer, //TODO: is this used anywhere?
		bool periodic
#ifdef QUICKSCHED
		, qsched *scheduler
#endif
		) : PseudoParticleContainer(orderOfExpansions),
			_leafContainer(nullptr),
			_wellSep(1) {
	_doNTLocal = true;
	_doNTGlobal = true;
	_periodicBC = periodic;
	_fuseGlobalCommunication = false;
	bool doDynamicAdjustment = true; //TODO: is this used anywhere?
#ifdef ENABLE_MPI
	int size;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//set avoidAllReduce depending on number of MPI ranks
	if(size >= 512){
		_avoidAllReduce = true;
	}
	else{
		if(_doNTGlobal and size >= 64){
			_avoidAllReduce = true;
		}
		else{
			_avoidAllReduce = false;
			if(size == 1){
				_avoidAllReduce = false;
			}
		}
	}
	if(size <= 64){ // in this case there might be doubling effects if not switched off (+y and -y neighbor the same -> twice added the same value in tower backcommunication for 1st local level)
		_doNTLocal = false;
	}
#else
	_avoidAllReduce = false;
#endif
	if(!_avoidAllReduce){
		_doNTGlobal = false;
	}
	//enable overlapping communication;
#if defined(ENABLE_MPI)
	_overlapComm = 1;
#else
	_overlapComm = 0;
#endif

#if defined(ENABLE_MPI)
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	if(!_overlapComm){
		std::vector<int> neigh = domainDecomp.getNeighbourRanks();
		for(int i = 0; i<6; ++i){
			_neighbours.push_back(neigh[i]);
		}
	}
	else{
		if(_avoidAllReduce){
			_allRanks = domainDecomp.getAllRanks();
		}
		std::vector<int> neigh = domainDecomp.getNeighbourRanksFullShell();
		for(int i = 0; i<26; ++i){
			_neighbours.push_back(neigh[i]);
			int myRank;
			MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
		}
		MPI_Barrier(MPI_COMM_WORLD);


	}
	_comm = domainDecomp.getCommunicator();
#endif
#if WIGNER == 1
	//global_log->error() << "not supported yet" << std::endl;
	mardyn_exit(-1);
#endif
#ifdef ENABLE_MPI
	/*
	//they are per default false
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_CELLS", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_GLOBAL", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_LOKAL", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_LOKAL", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_LOKAL", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_FAR_FIELD", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMMUNICATION_HALOS", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_HALO_GATHER", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_FMM_COMPLETE", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_CALCULATION", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_INIT", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_FINALIZE", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_EVAL_M", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_EVAL_LM", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE_ME", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_STOP_LEVEL", false);
	global_simulation->setSyncTimer("UNIFORM_PSEUDO_PARTICLE_CONTAINER_STOP_LEVEL", false);
	*/



#endif
	_leafContainer = new LeafNodesContainer(bBoxMin,
											bBoxMax,
											LJCellLength,
											LJSubdivisionFactor,
											periodic
#ifdef QUICKSCHED
											, scheduler
#endif
	);

	double cellLength[3];

	for (int i = 0; i < 3; i++) {
		cellLength[i] = _leafContainer->getCellLength()[i];
		_cellLength[i] = cellLength[i];
	}
#if defined(ENABLE_MPI)
	_bBoxMin[0] = bBoxMin[0];
	_bBoxMin[1] = bBoxMin[1];
	_bBoxMin[2] = bBoxMin[2];
#endif

	_globalNumCellsPerDim = domainLength[0] / cellLength[0];
	_maxLevel = log2(_globalNumCellsPerDim);
	mardyn_assert(_maxLevel == log2(domainLength[1] / cellLength[1]));
	mardyn_assert(_maxLevel == log2(domainLength[2] / cellLength[2]));
#if defined(ENABLE_MPI)
	int numProcessors;
	MPI_Comm_size(_comm,&numProcessors);
	_globalLevel = ceil(log2(numProcessors)/3.0);
	if(_globalLevel > _maxLevel){
		std::cout << "too many MPI ranks \n";
		mardyn_exit(-1);
	}
	//numProcessers has to be a power of 2
	mardyn_assert(pow(2,log2(numProcessors)) == numProcessors);
	//non overlap communication works only for powers of 8
	if(!_overlapComm){
		_numProcessorsPerDim = pow(2,log2(numProcessors)/3);
	}
	for(int i=0; i<3; i++){
		_numProcessorsPerDim[i] = rint(domainLength[i]/(bBoxMax[i] - bBoxMin[i]));
	}
	int maxProcsPerDim = std::max(_numProcessorsPerDim[0], std::max(_numProcessorsPerDim[1],_numProcessorsPerDim[2]));
	for(int i=0; i<3; i++){
		_numCellsOnGlobalLevel[i] = rint(maxProcsPerDim * 1.0 /_numProcessorsPerDim[i]);
	}
	if(_avoidAllReduce){
		if(_numCellsOnGlobalLevel[0] > 1 or _numCellsOnGlobalLevel[1] > 1 or _numCellsOnGlobalLevel[2] > 1){
			_importWholeGlobalRegion = true;
		}
		else{
			_importWholeGlobalRegion = false;
		}
	}
	else{
		_importWholeGlobalRegion = false;
	}
#else
	_globalLevel = _maxLevel;
#endif


	//allocate Multipole and Local particles
	int num_cells_in_level = 1;
	_mpCellGlobalTop.reserve(_globalLevel + 1);
#ifdef ENABLE_MPI
	_mpCellLocal.reserve(_maxLevel-_globalLevel);
#endif


	//num_cells_in_level_one_dim = 1;

	for (int n = 0; n <= _globalLevel; n++) {
		_mpCellGlobalTop.push_back(std::vector<MpCell>(num_cells_in_level, _maxOrd));
		num_cells_in_level *= 8;
	}
	//num_cells_in_level = 8;
	Vector3<int> num_cells_in_level_one_dim(_numCellsOnGlobalLevel * 2);
	//num_cells_in_level_one_dim = 2;

	int xHaloSize = 0;
	int yHaloSize = 0;
	int zHaloSize = 0;
	int yHaloSizeOverlap = 0;
	int zHaloSizeOverlap = 0;
	Vector3<int> edgeHaloSize = 0;
	int cornerHaloSize = 0;
	for (int n = _globalLevel + 1; n <= _maxLevel; n++) {
		_mpCellLocal.push_back(std::vector<MpCell>(static_cast<size_t>(num_cells_in_level_one_dim[0] + 4) * (num_cells_in_level_one_dim[1] + 4) * (num_cells_in_level_one_dim[2] + 4) , _maxOrd));
		xHaloSize += 2 * num_cells_in_level_one_dim[1] * num_cells_in_level_one_dim[2];
		yHaloSize += 2 * (num_cells_in_level_one_dim[0] + 4) * num_cells_in_level_one_dim[2];
		zHaloSize += 2 * (num_cells_in_level_one_dim[0] + 4) * (num_cells_in_level_one_dim[1] + 4);
		yHaloSizeOverlap += 2 * (num_cells_in_level_one_dim[0]) * num_cells_in_level_one_dim[2];
		zHaloSizeOverlap += 2 * (num_cells_in_level_one_dim[0]) * (num_cells_in_level_one_dim[1]);

		edgeHaloSize += num_cells_in_level_one_dim * 4;
		cornerHaloSize += 8;
		num_cells_in_level_one_dim *= 2;
	}

//	mardyn_assert(
//			num_cells_in_level
//					== _globalNumCellsPerDim * _globalNumCellsPerDim
//							* _globalNumCellsPerDim);


	// initalize centers and radii
	//num_cells_in_level = 1;
	Vector3<double> current_pos;
	Vector3<double> current_cell_length(domainLength);
	size_t num_cells_in_level_global = 1;
#if defined(ENABLE_MPI)
	Vector3<double> globalLevelCellLength = Vector3<double>(domainLength)
						* (1.0 / pow(2,_globalLevel));
	_processorPositionGlobalLevel[0] = rint(bBoxMin[0]/ globalLevelCellLength[0]);
	_processorPositionGlobalLevel[1] = rint(bBoxMin[1]/ globalLevelCellLength[1]);
	_processorPositionGlobalLevel[2] = rint(bBoxMin[2]/ globalLevelCellLength[2]);
#endif
	//initialization of global top part of tree
	for (int n = 0; n <= _globalLevel; ++n) {
		for (unsigned int z = 0; z < num_cells_in_level_global; ++z) {
			for (unsigned int y = 0; y < num_cells_in_level_global; ++y) {
				for (unsigned int x = 0; x < num_cells_in_level_global; ++x) {
					current_pos[0] = (x + 0.5) * current_cell_length[0];
					current_pos[1] = (y + 0.5) * current_cell_length[1];
					current_pos[2] = (z + 0.5) * current_cell_length[2];
					size_t cellIndex = ((z * num_cells_in_level_global + y)
							* num_cells_in_level_global) + x;
					_mpCellGlobalTop[n][cellIndex].multipole.setCenter(current_pos);
					_mpCellGlobalTop[n][cellIndex].multipole.setRadius(
							current_cell_length.L2Norm() * 0.5);

					_mpCellGlobalTop[n][cellIndex].local.setCenter(current_pos);
					_mpCellGlobalTop[n][cellIndex].local.setRadius(
							current_cell_length.L2Norm() * 0.5);
				}
			}
		}
		num_cells_in_level_global *= 2;
		//num_cells_in_level *= 8; //ToDo: is it needed anymore?
		current_cell_length = Vector3<double>(domainLength)
				* (1.0 / num_cells_in_level_global);
	}

	num_cells_in_level_one_dim = Vector3<int>(_numCellsOnGlobalLevel) * 2;
	//initialization of local tree
	for (int n = 0; n < _maxLevel-_globalLevel; ++n) {
		for (int z = -2; z < num_cells_in_level_one_dim[2]+2; ++z) {
			for (int y = -2; y < num_cells_in_level_one_dim[1]+2; ++y) {
				for (int x = -2; x < num_cells_in_level_one_dim[0]+2; ++x) {
					current_pos[0] = (x + 0.5) * current_cell_length[0] + bBoxMin[0];
					current_pos[1] = (y + 0.5) * current_cell_length[1] + bBoxMin[1];
					current_pos[2] = (z + 0.5) * current_cell_length[2] + bBoxMin[2];
					int cellIndex = (((z + 2) * (num_cells_in_level_one_dim[1] + 4) + y + 2)
							* (num_cells_in_level_one_dim[0] + 4)) + x + 2;
					_mpCellLocal[n][cellIndex].multipole.setCenter(current_pos);
					_mpCellLocal[n][cellIndex].multipole.setRadius(
							current_cell_length.L2Norm() * 0.5);

					_mpCellLocal[n][cellIndex].local.setCenter(current_pos);
					_mpCellLocal[n][cellIndex].local.setRadius(
							current_cell_length.L2Norm() * 0.5);
				}
			}
		}
		num_cells_in_level_one_dim *= 2;
		current_cell_length = current_cell_length * 0.5;
	}

	_domain = global_simulation->getDomain();
	_globalNumCells = pow(_globalNumCellsPerDim, 3);
#if defined(ENABLE_MPI)
	_globalLevelNumCells = pow(8,_globalLevel);
	_occVector = std::vector<int>(_globalLevelNumCells, 0);
#endif
	_coeffVectorLength = 0;
	_expansionSize = 0;
	for (int j = 0; j <= _maxOrd; j++) {
		for (int k = 0; k <= j; k++) {
			_expansionSize++;
		}
	}

#ifdef FMM_FFT
	//FFT objects initialization
	FFTSettings::autoSetting(_maxOrd); //auto set FFT acceleration's settings to optimal values
	//FFTSettings::TFMANAGER_VERBOSE = false;

	// TODO: move this some place sensible (Env var?, config?)
	taskModelTypesM2L taskModelM2L = CompleteTarget;
//    taskModelTypesM2L taskModelM2L = Pair2Way;
	if(taskModelM2L == Pair2Way)
		FFTSettings::USE_2WAY_M2L = true;
	else
		FFTSettings::USE_2WAY_M2L = false;

	FFTSettings::printCurrentOptions();
	if (FFTSettings::issetFFTAcceleration()) {
		_FFTAcceleration = FFTFactory::getFFTAccelerationAPI(_maxOrd);
		_FFT_TM = FFTFactory::getTransferFunctionManagerAPI(_maxOrd,
				_FFTAcceleration);
	} else {
		_FFTAcceleration = NULL;
		_FFT_TM = NULL;
	}
#endif  /* FMM_FFT */

#ifdef ENABLE_MPI
	//init neighborhoud coms and allreduce comms
	int coords[3];
	int myRank;
	MPI_Comm_rank(_comm,&myRank);
	MPI_Cart_coords(_comm, myRank, 3, coords);
	_neighbourhoodComms = new MPI_Comm[_globalLevel];
	_allReduceComms = new MPI_Comm[_globalLevel + 1];
	_allReduceComms[_globalLevel] = _comm;
	for(int i = _globalLevel-1; i>= 0; i--){
		int stride = pow(2,_globalLevel - i);
		MPI_Comm temp;
		int rowLength = pow(2,i);
		int colourDiv = (((coords[2] * _numCellsOnGlobalLevel[2]) /stride) * rowLength   + ((coords[1] * _numCellsOnGlobalLevel[1]) / stride)) * rowLength + ((coords[0] * _numCellsOnGlobalLevel[0]) / stride);
		int colourStride = (((coords[2] * _numCellsOnGlobalLevel[2]) % stride) * stride   + ((coords[1] * _numCellsOnGlobalLevel[1]) %  stride)) * stride + ((coords[0] * _numCellsOnGlobalLevel[0]) % stride);
		//calculate all processors that have to contribute in "global" allreduce after all levels until stop level i have been calculated
		MPI_Comm_split(_comm, colourStride, 0, &_allReduceComms[i]);
		//calculate neighbourhood cooms for avoiding allReduce
		MPI_Comm_split(_comm, colourDiv, 0, &temp);
		stride /= 2;
		colourStride = (((coords[2] * _numCellsOnGlobalLevel[2]) % stride) * stride   + ((coords[1] * _numCellsOnGlobalLevel[1]) %  stride)) * stride + ((coords[0] * _numCellsOnGlobalLevel[0]) % stride);
		MPI_Comm_split(temp, colourStride, 0, &_neighbourhoodComms[i]);
		int size2;
		MPI_Comm_size(_neighbourhoodComms[i], &size2);
		if(size2 > 8){ //neighbourhood comms need to have size 8
			std::cout << "Error wrong communicator \n";
			mardyn_exit(1);
		}
	}
#endif

#if defined(ENABLE_MPI)
	xHaloSize *= _expansionSize;
	yHaloSize *= _expansionSize;
	zHaloSize *= _expansionSize;
	yHaloSizeOverlap *= _expansionSize;
	zHaloSizeOverlap *= _expansionSize;
	cornerHaloSize *= _expansionSize;
	edgeHaloSize *= _expansionSize;
	Vector3<int> areaHaloSize;
	areaHaloSize[0] = xHaloSize;
	areaHaloSize[1] = yHaloSizeOverlap;
	areaHaloSize[2] = zHaloSizeOverlap;
	if(!_overlapComm){
		_multipoleBuffer = new HaloBufferNoOverlap<double>(xHaloSize * 2,yHaloSize * 2,zHaloSize * 2);
		_multipoleRecBuffer = new HaloBufferNoOverlap<double>(xHaloSize * 2,yHaloSize * 2,zHaloSize * 2);
	}
	else{
		std::vector<int> areaNeighbours;
		std::vector<int> edgesNeighbours;
		std::vector<int> cornerNeighbours;

		for(int i = 0; i < 6 ; i++){
			areaNeighbours.push_back(_neighbours[i]);
		}
		for(int i = 6; i < 18 ; i++){
			edgesNeighbours.push_back(_neighbours[i]);
		}
		for(int i = 18; i < 26 ; i++){
			cornerNeighbours.push_back(_neighbours[i]);
		}
		_multipoleBufferOverlap = new HaloBufferOverlap<double>(areaHaloSize * 2,edgeHaloSize * 2, cornerHaloSize * 2, _comm, areaNeighbours, edgesNeighbours, cornerNeighbours, 1, _doNTLocal);
		_multipoleRecBufferOverlap = new HaloBufferOverlap<double>(areaHaloSize * 2,edgeHaloSize * 2, cornerHaloSize * 2, _comm, areaNeighbours, edgesNeighbours, cornerNeighbours, 0, _doNTLocal);
		if(_avoidAllReduce){ //initialize global buffers and optimize stopping level
			int importVolumePerLevel;
			if(_fuseGlobalCommunication){
				importVolumePerLevel = 26;
			}
			else{
				importVolumePerLevel = _importWholeGlobalRegion? 216 : 189;
			}
			//import on every level
			int importVolume = importVolumePerLevel * _globalLevel;
			if(_fuseGlobalCommunication) importVolume -= importVolumePerLevel; //at level 1 there is no need for communication due to fuse step
			int cellsPerMessage = (_fuseGlobalCommunication)? 8 : 1;
			int numberOfGlobalLevel = _globalLevel;
			if(_fuseGlobalCommunication) numberOfGlobalLevel--;
			_multipoleBufferOverlapGlobal =  new HaloBufferOverlap<double>(Vector3<int>(cellsPerMessage *_expansionSize * 2), Vector3<int>(0), 0, _comm, areaNeighbours, edgesNeighbours, cornerNeighbours, 1, _doNTGlobal, importVolume, 0, 0, _allRanks,_numCellsOnGlobalLevel,_fuseGlobalCommunication);
			_multipoleRecBufferOverlapGlobal = new HaloBufferOverlap<double>(Vector3<int>(cellsPerMessage *_expansionSize * 2), Vector3<int>(0), 0, _comm, areaNeighbours, edgesNeighbours, cornerNeighbours, 0, _doNTGlobal, importVolume, 0, 0, _allRanks,_numCellsOnGlobalLevel,_fuseGlobalCommunication);
			_multipoleBufferOverlapGlobal->setNumberOfGlobalLevelsInBuffer(numberOfGlobalLevel);
			_multipoleRecBufferOverlapGlobal->setNumberOfGlobalLevelsInBuffer(numberOfGlobalLevel);

			//get optimal stopping level for change from local allreduce to global allreduce
			if(doDynamicAdjustment){
				_stopLevel = optimizeAllReduce(/*ljContainer*/);
			}
			else{
				_stopLevel = 1;
			}
			if(myRank == 0){
				std::cout << "optimal stop level = " << _stopLevel << "\n";
			}
			delete _multipoleBufferOverlapGlobal;
			delete _multipoleRecBufferOverlapGlobal;
			if(_stopLevel <= _globalLevel){ //reinitialize buffers with optimized stoplevel
				importVolume = importVolumePerLevel * (_globalLevel - _stopLevel + 1);
				if(_fuseGlobalCommunication and _stopLevel == 1) importVolume -= importVolumePerLevel;
				_multipoleBufferOverlapGlobal =  new HaloBufferOverlap<double>(Vector3<int>(cellsPerMessage * _expansionSize * 2), Vector3<int>(0), 0, _comm, areaNeighbours, edgesNeighbours, cornerNeighbours, 1, _doNTGlobal, importVolume, 0, 0, _allRanks,_numCellsOnGlobalLevel,_fuseGlobalCommunication);
				_multipoleRecBufferOverlapGlobal = new HaloBufferOverlap<double>(Vector3<int>(cellsPerMessage * _expansionSize * 2), Vector3<int>(0), 0, _comm, areaNeighbours, edgesNeighbours, cornerNeighbours, 0, _doNTGlobal, importVolume, 0, 0, _allRanks,_numCellsOnGlobalLevel,_fuseGlobalCommunication);
				numberOfGlobalLevel = _globalLevel + 1 - _stopLevel;
				if(_fuseGlobalCommunication  and _stopLevel == 1) numberOfGlobalLevel--;
				_multipoleBufferOverlapGlobal->setNumberOfGlobalLevelsInBuffer(numberOfGlobalLevel);
				_multipoleRecBufferOverlapGlobal->setNumberOfGlobalLevelsInBuffer(numberOfGlobalLevel);

				if(_globalLevel >= 1 and not (_globalLevel == 1 and _fuseGlobalCommunication)){
					_multipoleRecBufferOverlapGlobal->communicateGlobalLevels(_globalLevel, _stopLevel);
				}
			}
			else{
				_avoidAllReduce = false;
				_doNTGlobal = false;
				_fuseGlobalCommunication = false;
			}
		}
		_multipoleRecBufferOverlap->communicate(false);

		//take care that every process has initiated receive in first iteration
		MPI_Barrier(_comm);
	}
#endif

#if defined(ENABLE_MPI)
	int numCells = 0;
	int startLevel;
	if(_avoidAllReduce){
		startLevel = _stopLevel - 1;
	}
	else{
		startLevel = _globalLevel;
	}
	for(int i = startLevel; i >=0 ; i--){
		numCells += pow(8,i);
	}
	_coeffVectorLength = _expansionSize*numCells;
#endif

	_coeffVector = std::vector<double>(_coeffVectorLength * 2, 0.0);
	Log::global_log->info() << "UniformPseudoParticleContainer: coeffVectorLength="
							<< _coeffVectorLength << " Size of MPI Buffers is "
							<< (8 * (_coeffVectorLength * 2 + _globalNumCells)
									/ (1024.0 * 1024.0)) << " MB;" << std::endl;
	Log::global_log->info() << "UniformPseudoParticleContainer: globalLevel = "
							<< _globalLevel << " maxLevel= " <<_maxLevel << std::endl;

#ifdef QUICKSCHED
	// Quicksched Task and resource generation generation
	generateResources(scheduler);
	generateL2PTasks(scheduler);
	generateP2PTasks(scheduler);
	generateL2LTasks(scheduler);
	generateM2LTasks(scheduler, taskModelM2L);
	generateP2MTasks(scheduler);
	generateM2MTasks(scheduler);
#endif
//reset timers
#ifdef ENABLE_MPI
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_CELLS");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE_ME");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_GLOBAL");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_LOKAL");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_LOKAL");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_LOKAL");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROCESS_FAR_FIELD");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMMUNICATION_HALOS");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_HALO_GATHER");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_FMM_COMPLETE");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_CALCULATION");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_INIT");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_FINALIZE");
	global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL");
#endif

}

UniformPseudoParticleContainer::~UniformPseudoParticleContainer() {
	delete _leafContainer;
#if defined(ENABLE_MPI)
	if(!_overlapComm){
		delete _multipoleBuffer;
		delete _multipoleRecBuffer;
	}
	else{
		delete _multipoleBufferOverlap;
		delete _multipoleRecBufferOverlap;
		if(_avoidAllReduce){
			delete _multipoleBufferOverlapGlobal;
			delete _multipoleRecBufferOverlapGlobal;
		}
	}
#endif

#ifdef FMM_FFT
	if (FFTSettings::issetFFTAcceleration()) {
		delete _FFT_TM;
		delete _FFTAcceleration;
	}
#endif  /* FMM_FFT */
}

#ifdef QUICKSCHED
void UniformPseudoParticleContainer::generateResources(qsched *scheduler) {

	// resources for multipoles
	for (int level = 0; level <= _maxLevel; ++level) {
		for (unsigned int multipoleId = 0; multipoleId < _mpCellGlobalTop[level].size(); ++multipoleId) {
			_mpCellGlobalTop[level][multipoleId]._resIdLocal     = qsched_addres(scheduler,
																				 qsched_owner_none,
																				 qsched_res_none);
			_mpCellGlobalTop[level][multipoleId]._resIdMultipole = qsched_addres(scheduler,
																				 qsched_owner_none,
																				 qsched_res_none);
		}
	}
	// resources for lj cells
	for (unsigned int i = 0; i < _leafContainer->getCells().size(); ++i) {
		_leafContainer->getCells()[i].setResourceId(
				qsched_addres(scheduler, qsched_owner_none, qsched_res_none)
		);
	}
}

void UniformPseudoParticleContainer::generateP2MTasks(qsched *scheduler) {
	struct qsched_payload payload;
	payload.uniformPseudoParticleContainer = this;

	long sourceId,
		 targetId;
	int  idP2M,
		 targetCoords[3],
		 cellsPerDim = 1 << _maxLevel;     // in cells

	for (targetCoords[0] = 0; targetCoords[0] < cellsPerDim; ++targetCoords[0]) {
		for (targetCoords[1] = 0; targetCoords[1] < cellsPerDim; ++targetCoords[1]) {
			for (targetCoords[2] = 0; targetCoords[2] < cellsPerDim; ++targetCoords[2]) {
				sourceId = _leafContainer->cellIndexOf3DIndex(targetCoords[0]+1,
															  targetCoords[1]+1,
															  targetCoords[2]+1);
				payload.currentMultipole = sourceId;
				idP2M = qsched_addtask(scheduler,
									   FastMultipoleMethod::P2MCompleteCell,
									   qsched_flag_none,
									   &payload,
									   sizeof(payload),
									   1);

				targetId = (targetCoords[2] * cellsPerDim + targetCoords[1])
						   * cellsPerDim
						   + targetCoords[0];
				_mpCellGlobalTop[_maxLevel][targetId]._taskIdP2M = idP2M;
				qsched_addunlock(scheduler,
								 idP2M,
								 _mpCellGlobalTop[_maxLevel][targetId]._taskIdM2LInit
				);
			}
		}
	}
}

void UniformPseudoParticleContainer::generateM2MTasks(qsched *scheduler) {
	struct qsched_payload payload;
	payload.uniformPseudoParticleContainer = this;

	int idM2M,
		sourceId,
		cellsPerDim  = 1 << (_maxLevel - 1),       // in cells
		cellsInLevel = 1 << (3 * (_maxLevel - 1)), // == cellsPerDim^3
		sourceCoords[3],
		targetCoords[3];

	for (int level = _maxLevel - 1; level > 0; --level) {
		payload.currentLevel = level;
		payload.currentEdgeLength = cellsPerDim;
		for (int targetId = 0; targetId < cellsInLevel; ++targetId) {
			payload.currentMultipole = targetId;

			idM2M = qsched_addtask(scheduler,
								   FastMultipoleMethod::M2MCompleteCell,
								   task_flag_none,
								   &payload,
								   sizeof(payload),
								   1);

			_mpCellGlobalTop[level][targetId]._taskIdM2M = idM2M;

			qsched_addunlock(scheduler,
							 idM2M,
							 _mpCellGlobalTop[level][targetId]._taskIdM2LInit);

			int cellsPerDim2 = cellsPerDim * cellsPerDim;
			targetCoords[2] = targetId / cellsPerDim2;
			targetCoords[1] = (targetId - (targetCoords[2] * cellsPerDim2)) / cellsPerDim;
			targetCoords[0] = targetId - (targetCoords[2] * cellsPerDim2) - (targetCoords[1] * cellsPerDim);

			for (int i = 0; i < 8; ++i) {
				sourceCoords[0] = targetCoords[0] * 2;
				sourceCoords[1] = targetCoords[1] * 2;
				sourceCoords[2] = targetCoords[2] * 2;

				if (IsOdd(i    )) ++sourceCoords[0];
				if (IsOdd(i / 2)) ++sourceCoords[1];
				if (IsOdd(i / 4)) ++sourceCoords[2];

				sourceId = (sourceCoords[2] * cellsPerDim * 2 + sourceCoords[1]) * cellsPerDim * 2
						   + sourceCoords[0];

				if (level != _maxLevel - 1) {
					qsched_addunlock(scheduler,
									 _mpCellGlobalTop[level + 1][sourceId]._taskIdM2M,
									 _mpCellGlobalTop[level][targetId]._taskIdM2M);
				} else { // level == _maxLevel - 1
					qsched_addunlock(scheduler,
									 _mpCellGlobalTop[level + 1][sourceId]._taskIdP2M,
									 _mpCellGlobalTop[level][targetId]._taskIdM2M);
				}
			}
		}
		cellsPerDim >>= 1;
		cellsInLevel >>= 3;
	}
}

void UniformPseudoParticleContainer::generateM2LTasks(qsched *scheduler, taskModelTypesM2L taskModelM2L) {

	qsched_task_t         idInit,
						  idFin,
						  idM2L;
	struct qsched_payload payload;
	payload.uniformPseudoParticleContainer = this;

	switch (taskModelM2L) {
		case CompleteTarget: {
			int currentEdgeLength = 2,
				cellsInLevel      = 8;
			for (int currentLevel = 1; currentLevel <= _maxLevel; ++currentLevel) {

				payload.currentLevel      = currentLevel;
				payload.currentEdgeLength = currentEdgeLength;

				// generate tasks
				for (int multipoleId = 0; multipoleId < cellsInLevel; ++multipoleId) {
					payload.currentMultipole = multipoleId;
					idInit = qsched_addtask(scheduler,
											FastMultipoleMethod::M2LInitializeSource,
											task_flag_none,
											&payload,
											sizeof(payload),
											1);
					_mpCellGlobalTop[currentLevel][multipoleId]._taskIdM2LInit = idInit;
					idM2L = qsched_addtask(scheduler,
										   FastMultipoleMethod::M2LTranslation,
										   task_flag_none,
										   &payload,
										   sizeof(payload),
										   1);
					_mpCellGlobalTop[currentLevel][multipoleId]._taskIdM2LCalc = idM2L;
					// Translation tasks includes finalize
					_mpCellGlobalTop[currentLevel][multipoleId]._taskIdM2LFin  = idM2L;
				}

				// generate dependencies
				// TODO: think of smarter iteration scheme to merge this and task loop
				for (int targetId = 0; targetId < cellsInLevel; ++targetId) {
					int sourceCoords[3],
						targetCoords[3];
					int currentEdgeLength2 = currentEdgeLength * currentEdgeLength;
					targetCoords[2] = targetId / currentEdgeLength2;
					targetCoords[1] = (targetId - (targetCoords[2] * currentEdgeLength2)) / currentEdgeLength;
					targetCoords[0] = targetId - (targetCoords[2] * currentEdgeLength2) - (targetCoords[1] * currentEdgeLength);

					// cell from where sources will be calculated
					int      sourceRootCoords[3];
					// find start corner for stencil
					for (int i = 0; i < 3; ++i) {
						if (IsOdd(targetCoords[i]))
							sourceRootCoords[i] = targetCoords[i] - 3;
						else
							sourceRootCoords[i] = targetCoords[i] - 2;
					}

					// finalizing (here done in Calc) unlocks L2L. No L2Ls on maxLevel.
					if(currentLevel < _maxLevel){
						qsched_addunlock(scheduler,
										 _mpCellGlobalTop[currentLevel][targetId]._taskIdM2LFin,
										 _mpCellGlobalTop[currentLevel][targetId]._taskIdL2L);
					} else { // currentLevel == _maxLevel
						qsched_addunlock(scheduler,
										 _mpCellGlobalTop[currentLevel][targetId]._taskIdM2LFin,
										 _mpCellGlobalTop[currentLevel][targetId]._taskIdL2P);
					}
					qsched_addlock(scheduler,
								   _mpCellGlobalTop[currentLevel][targetId]._taskIdM2LCalc,
								   _mpCellGlobalTop[currentLevel][targetId]._resIdLocal);

					// iterate through 6x6x6 grid around target cell with periodic boundaries
					for (int x = 0; x < 6; ++x) {
						sourceCoords[0] = sourceRootCoords[0] + x;
						if (sourceCoords[0] < 0)
							sourceCoords[0] += currentEdgeLength;
						else if (sourceCoords[0] >= currentEdgeLength)
							sourceCoords[0] -= currentEdgeLength;
						for (int y = 0; y < 6; ++y) {
							sourceCoords[1] = sourceRootCoords[1] + y;
							if (sourceCoords[1] < 0)
								sourceCoords[1] += currentEdgeLength;
							else if (sourceCoords[1] >= currentEdgeLength)
								sourceCoords[1] -= currentEdgeLength;
							for (int z = 0; z < 6; ++z) {
								sourceCoords[2] = sourceRootCoords[2] + z;
								if (sourceCoords[2] < 0)
									sourceCoords[2] += currentEdgeLength;
								else if (sourceCoords[2] >= currentEdgeLength)
									sourceCoords[2] -= currentEdgeLength;
								// skip cells that are too near
								if (std::abs(sourceRootCoords[0] + x - targetCoords[0]) < 2
									&& std::abs(sourceRootCoords[1] + y - targetCoords[1]) < 2
									&& std::abs(sourceRootCoords[2] + z - targetCoords[2]) < 2)
									continue;

								int sourceId = (sourceCoords[2] * currentEdgeLength + sourceCoords[1])
											   * currentEdgeLength
											   + sourceCoords[0];

								qsched_addunlock(scheduler,
												 _mpCellGlobalTop[currentLevel][sourceId]._taskIdM2LInit,
												 _mpCellGlobalTop[currentLevel][targetId]._taskIdM2LCalc);
							}
						}
					}
				}

				currentEdgeLength <<= 1;
				// cellsInLevel = currentEdgeLength * currentEdgeLength * currentEdgeLength
				cellsInLevel <<= 3;

			}
			break;
		} /* CompleteTarget */
		case Pair2Way : {
			// TODO: might be broken since downward pass was added
			int currentEdgeLength = 2,
				cellsInLevel      = 8;
			for (int currentLevel = 1; currentLevel <= _maxLevel; ++currentLevel) {

				payload.currentLevel      = currentLevel;
				payload.currentEdgeLength = currentEdgeLength;

				// generate resource ids, init and finalize tasks
				for (int multipoleId = 0; multipoleId < cellsInLevel; ++multipoleId) {
					payload.currentMultipole = multipoleId;
					idInit = qsched_addtask(scheduler,
											FastMultipoleMethod::M2LInitializeCell,
											task_flag_none,
											&payload,
											sizeof(payload),
											1);
					_mpCellGlobalTop[currentLevel][multipoleId]._taskIdM2LInit = idInit;
					idFin = qsched_addtask(scheduler,
										   FastMultipoleMethod::M2LFinalizeCell,
										   task_flag_none,
										   &payload,
										   sizeof(payload),
										   1);
					_mpCellGlobalTop[currentLevel][multipoleId]._taskIdM2LFin = idFin;
					// finalizing unlocks L2L
					qsched_addunlock(scheduler,
									 _mpCellGlobalTop[currentLevel][multipoleId]._taskIdM2LFin,
									 _mpCellGlobalTop[currentLevel][multipoleId]._taskIdL2L);
				}
				// generate m2l tasks and dependencies
				for (int targetId = 0; targetId < cellsInLevel; ++targetId) {
					payload.currentMultipole = targetId;
					int sourceCoords[3],
						targetCoords[3];
					int currentEdgeLength2 = currentEdgeLength * currentEdgeLength;
					targetCoords[0] = targetId / currentEdgeLength2;
					targetCoords[1] = (targetId - (targetCoords[0] * currentEdgeLength2)) / currentEdgeLength;
					targetCoords[2] = targetId - (targetCoords[0] * currentEdgeLength2) - (targetCoords[1] * currentEdgeLength);
					// cell from where sources will be calculated
					int sourceRootCoords[3];

					// find start corner for stencil
					for (int i = 0; i < 3; ++i) {
						if (IsOdd(targetCoords[i]))
							sourceRootCoords[i] = targetCoords[i] - 3;
						else
							sourceRootCoords[i] = targetCoords[i] - 2;
					}

					// iterate through 6x6x6 grid around target cell with periodic boundaries
					for (int x = 0; x < 6; ++x) {
						sourceCoords[0] = sourceRootCoords[0] + x;
						if (sourceCoords[0] < 0)
							sourceCoords[0] += currentEdgeLength;
						else if (sourceCoords[0] >= currentEdgeLength)
							sourceCoords[0] -= currentEdgeLength;
						for (int y = 0; y < 6; ++y) {
							sourceCoords[1] = sourceRootCoords[1] + y;
							if (sourceCoords[1] < 0)
								sourceCoords[1] += currentEdgeLength;
							else if (sourceCoords[1] >= currentEdgeLength)
								sourceCoords[1] -= currentEdgeLength;
							for (int z = 0; z < 6; ++z) {
								sourceCoords[2] = sourceRootCoords[2] + z;
								if (sourceCoords[2] < 0)
									sourceCoords[2] += currentEdgeLength;
								else if (sourceCoords[2] >= currentEdgeLength)
									sourceCoords[2] -= currentEdgeLength;
								// skip cells that are too near
								if (std::abs(sourceRootCoords[0] + x- targetCoords[0]) < 2
									&& std::abs(sourceRootCoords[1] + y - targetCoords[1]) < 2
									&& std::abs(sourceRootCoords[2] + z - targetCoords[2]) < 2)
									continue;

								int sourceId = (sourceCoords[2] * currentEdgeLength + sourceCoords[1])
											   * currentEdgeLength
											   + sourceCoords[0];
								// FIXME: this prevents double calculation but should be replaced by smarter domain iteration
								if (sourceId > targetId)
									continue;
								payload.sourceMultipole = sourceId;

								idM2L = qsched_addtask(scheduler,
													   FastMultipoleMethod::M2LPair2Way,
													   task_flag_none,
													   &payload,
													   sizeof(payload),
													   1);

								// every calculation task depends on source and target initialization
								qsched_addunlock(scheduler,
												 _mpCellGlobalTop[currentLevel][targetId]._taskIdM2LInit,
												 idM2L);
								qsched_addunlock(scheduler,
												 _mpCellGlobalTop[currentLevel][sourceId]._taskIdM2LInit,
												 idM2L);
								// every calculation tasks brings one one step closer to finalizing target and source
								qsched_addunlock(scheduler,
												 idM2L,
												 _mpCellGlobalTop[currentLevel][targetId]._taskIdM2LFin);
								qsched_addunlock(scheduler,
												 idM2L,
												 _mpCellGlobalTop[currentLevel][sourceId]._taskIdM2LFin);
								// every calculation needs to lock source and target
								qsched_addlock(scheduler,
											   idM2L,
											   _mpCellGlobalTop[currentLevel][targetId]._resIdLocal);
								qsched_addlock(scheduler,
											   idM2L,
											   _mpCellGlobalTop[currentLevel][sourceId]._resIdLocal);
							}
						}
					}
				}

				currentEdgeLength <<= 1;
				// cellsInLevel = currentEdgeLength * currentEdgeLength * currentEdgeLength
				cellsInLevel <<= 3;
			}
			break;
		} /* Pair2Way */
	}
}

void UniformPseudoParticleContainer::generateL2LTasks(qsched *scheduler) {
	struct qsched_payload payload;
	payload.uniformPseudoParticleContainer = this;

	int idL2L,
		targetId,
		cellsPerDim  = 1 << (_maxLevel - 1),
		cellsInLevel = 1 << (3 * (_maxLevel - 1)), // == cellsPerDim^3
		sourceCoords[3],
		targetCoords[3];

	for (int level = _maxLevel-1; level > 0; --level) {

		payload.currentLevel = level;
		payload.currentEdgeLength = cellsPerDim;

		for (int sourceId = 0; sourceId < cellsInLevel; ++sourceId) {

			payload.sourceMultipole = sourceId;

			idL2L = qsched_addtask(scheduler,
								   FastMultipoleMethod::L2LCompleteCell,
								   task_flag_none,
								   &payload,
								   sizeof(payload),
								   1);

			_mpCellGlobalTop[level][sourceId]._taskIdL2L = idL2L;

			int cellsPerDim2 = cellsPerDim * cellsPerDim;
			sourceCoords[2] = sourceId / cellsPerDim2;
			sourceCoords[1] = (sourceId - (sourceCoords[2] * cellsPerDim2)) / cellsPerDim;
			sourceCoords[0] = sourceId - (sourceCoords[2] * cellsPerDim2) - (sourceCoords[1] * cellsPerDim);


			for (int i = 0; i < 8; ++i) {
				targetCoords[2] = 2 * sourceCoords[2];
				targetCoords[1] = 2 * sourceCoords[1];
				targetCoords[0] = 2 * sourceCoords[0];

				if (IsOdd(i    )) targetCoords[0] = targetCoords[0] + 1;
				if (IsOdd(i / 2)) targetCoords[1] = targetCoords[1] + 1;
				if (IsOdd(i / 4)) targetCoords[2] = targetCoords[2] + 1;

				targetId = (targetCoords[2] * cellsPerDim * 2 + targetCoords[1]) * cellsPerDim * 2
						   + targetCoords[0];

				// prevent concurrent M2L to this target
				qsched_addlock(scheduler,
							   _mpCellGlobalTop[level][sourceId]._taskIdL2L,
							   _mpCellGlobalTop[level + 1][targetId]._resIdLocal);
				if ((level + 1) != _maxLevel) {
					qsched_addunlock(scheduler,
									 _mpCellGlobalTop[level][sourceId]._taskIdL2L,
									 _mpCellGlobalTop[level + 1][targetId]._taskIdL2L);
				} else {
					qsched_addunlock(scheduler,
									 _mpCellGlobalTop[level][sourceId]._taskIdL2L,
									 _mpCellGlobalTop[level + 1][targetId]._taskIdL2P);
				}
			}
		}

		// cellsPerDim / 2;
		cellsPerDim >>= 1;
		// cellsInLevel / 8;
		cellsInLevel >>= 3;
	}
}

void UniformPseudoParticleContainer::generateL2PTasks(qsched *scheduler) {
	long sourceId,
		 targetId;
	int  idL2P,
		 sourceCoords[3],
		 cellsPerDim = 1 << _maxLevel;     // in cells

	struct qsched_payload payload;

	payload.uniformPseudoParticleContainer = this;

	for (sourceCoords[0] = 0; sourceCoords[0] < cellsPerDim; ++sourceCoords[0]) {
		for (sourceCoords[1] = 0; sourceCoords[1] < cellsPerDim; ++sourceCoords[1]) {
			for (sourceCoords[2] = 0; sourceCoords[2] < cellsPerDim; ++sourceCoords[2]) {
				targetId = _leafContainer->cellIndexOf3DIndex(sourceCoords[0]+1,
															  sourceCoords[1]+1,
															  sourceCoords[2]+1);
				payload.currentMultipole = targetId;
				idL2P = qsched_addtask(scheduler,
									   FastMultipoleMethod::L2PCompleteCell,
									   qsched_flag_none,
									   &payload,
									   sizeof(payload),
									   1);
				qsched_addlock(scheduler,
							   idL2P,
							   _leafContainer->getCells()[targetId].getTaskData()._resourceId);
				sourceId = (sourceCoords[2] * cellsPerDim + sourceCoords[1])
						   * cellsPerDim
						   + sourceCoords[0];
				_mpCellGlobalTop[_maxLevel][sourceId]._taskIdL2P = idL2P;
			}
		}
	}
}

void UniformPseudoParticleContainer::generateP2PTasks(qsched *scheduler) {
	long cellIndex;
	int  multipolesPerDim = 1 << _maxLevel,
		 cellsPerDim      = multipolesPerDim + 2;
	//TODO: Make taskBlockSize dynamic / autotuning / xml
	struct qsched_payload payload;
	payload.taskBlockSize[0] = 2;
	payload.taskBlockSize[1] = 2;
	payload.taskBlockSize[2] = 2;
	payload.leafNodesContainer = _leafContainer;


	Log::global_log->info() << "Generating P2P task ids" << std::endl;
	for (auto z = 0; z < cellsPerDim; ++z) {
		for (auto y = 0; y < cellsPerDim; ++y) {
			for (auto x = 0; x < cellsPerDim; ++x) {
				cellIndex = _leafContainer->cellIndexOf3DIndex(x, y, z);
				auto *cell = &_leafContainer->getCells()[cellIndex];
				// only create tasks with offset blocksize-1.
				// -1 because they need to overlap
				// skip tasks for rear halo layer as they would only contain halo cells
				if ((z % (payload.taskBlockSize[2] - 1) == 0
					 && y % (payload.taskBlockSize[1] - 1) == 0
					 && x % (payload.taskBlockSize[0] - 1) == 0)
					&&
					(x < cellsPerDim - 1
					 && y < cellsPerDim - 1
					 && z < cellsPerDim - 1)) {
					// P2P TASK
					// also save the pointers as long
					payload.cell.coordinates[0] = x;
					payload.cell.coordinates[1] = y;
					payload.cell.coordinates[2] = z;
					cell->setP2PId(qsched_addtask(scheduler,
												 FastMultipoleMethod::P2Pc08StepBlock,
												 task_flag_none,
												 &payload,
												 sizeof(payload),
												 1)
					);
				}
				// Pre and post tasks for every cell
				// PREPROCESS TASK
				payload.cell.pointer = cell;
				cell->setPreprocessId(qsched_addtask(scheduler,
													FastMultipoleMethod::P2PPreprocessSingleCell,
													task_flag_none,
													&payload,
													sizeof(payload),
													1)
				);
				// POSTPROCESS TASK
				cell->setPostprocessId(qsched_addtask(scheduler,
													 FastMultipoleMethod::P2PPostprocessSingleCell,
													 task_flag_none,
													 &payload,
													 sizeof(payload),
													 1)
				);
			} /* end for-x */
		} /* end for-y*/
	} /* end for-z */

	// set dependencies
	Log::global_log->info() << "Setting P2P task dependencies" << std::endl;
	for (auto z = 0; z < cellsPerDim - 1; z += payload.taskBlockSize[2] - 1) {
		for (auto y = 0; y < cellsPerDim - 1; y += payload.taskBlockSize[1] - 1) {
			for (auto x = 0; x < cellsPerDim - 1; x += payload.taskBlockSize[0] - 1) {
				cellIndex = _leafContainer->cellIndexOf3DIndex(x, y, z);
				auto *baseCell = &_leafContainer->getCells()[cellIndex];

				for (auto k = 0; k < payload.taskBlockSize[2]
								 && z + k < cellsPerDim; ++k) {
					for (auto j = 0; j < payload.taskBlockSize[1]
									 && y + j < cellsPerDim; ++j) {
						for (auto i = 0; i < payload.taskBlockSize[0]
										 && x + i < cellsPerDim; ++i) {
							cellIndex = _leafContainer->cellIndexOf3DIndex(x + k, y + j, z + i);
							auto *targetCell = &_leafContainer->getCells()[cellIndex];
							// create locks for only for resources at edges
							if(i == payload.taskBlockSize[0] - 1
							   || j == payload.taskBlockSize[1] - 1
							   || k == payload.taskBlockSize[2] - 1){
								qsched_addlock(scheduler,
											   baseCell->getTaskData()._P2PId,
											   targetCell->getTaskData()._resourceId);
							}
							// 8 preprocess (partly) unlock one P2P
							qsched_addunlock(scheduler,
											 targetCell->getTaskData()._preprocessId,
											 baseCell->getTaskData()._P2PId);
							// every P2P (partly) unlocks 8 postprocess
							qsched_addunlock(scheduler,
											 baseCell->getTaskData()._P2PId,
											 targetCell->getTaskData()._postprocessId);
						}
					}
				}
			} /* end for-x */
		} /* end for-y*/
	} /* end for-z */
}
#endif

void UniformPseudoParticleContainer::build(ParticleContainer* pc) {
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_FMM_COMPLETE");
	_leafContainer->clearParticles();
	for(auto tM = pc->iterator(ParticleIterator::ALL_CELLS); tM.isValid(); ++tM) {
		_leafContainer->addParticle(*tM);
	}
}

int UniformPseudoParticleContainer::optimizeAllReduce(/*ParticleContainer* ljContainer*/){
#if defined(ENABLE_MPI)
//	ljContainer->deleteOuterParticles();
//	ljContainer->update();
//	ljContainer->updateMoleculeCaches();
	P2MCellProcessor * _P2MProcessor = new P2MCellProcessor(this);
	L2PCellProcessor * _L2PProcessor = new L2PCellProcessor(this);
	VectorizedChargeP2PCellProcessor *_P2PProcessor = new VectorizedChargeP2PCellProcessor(
				*(global_simulation->getDomain()));
	double minTime = pow(2,100);
	int bestStopLevel = 1;
	std::vector<double> timeArray(_globalLevel);
	for(int stopLevel = 1; stopLevel <= _globalLevel + 1; stopLevel++){ //iterate over possible stopping level
		_stopLevel = stopLevel;
		if(_globalLevel >= 1 and not(_globalLevel == 1 and _fuseGlobalCommunication)){
			_multipoleRecBufferOverlapGlobal->communicateGlobalLevels(_globalLevel, _stopLevel);
		}
		_multipoleRecBufferOverlap->communicate(false);

		MPI_Barrier(_comm);
		int numCells = 0;
		int startLevel;
		if(_avoidAllReduce){
			startLevel = _stopLevel - 1;
		}
		else{
			startLevel = _globalLevel;
		}
		for(int i = startLevel; i >=0 ; i--){
			numCells += pow(8,i);
		}
		_coeffVectorLength = _expansionSize*numCells;

		_coeffVector = std::vector<double>(_coeffVectorLength * 2);

		_leafContainer->clearParticles();

		/*Molecule* tM;
		for(tM = ljContainer->begin(); tM != ljContainer->end(); tM = ljContainer->next()) {
			_leafContainer->addParticle(*tM);
		}*/
//		 clear expansions
		clear();

		global_simulation->timers()->reset("UNIFORM_PSEUDO_PARTICLE_CONTAINER_STOP_LEVEL");
		global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_STOP_LEVEL");

		// P2M, M2P
		upwardPass(_P2MProcessor);

		horizontalPass(_P2PProcessor);
		global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_STOP_LEVEL");
		timeArray[stopLevel] = global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_STOP_LEVEL");
		MPI_Allreduce(MPI_IN_PLACE,&timeArray[stopLevel],1,MPI_DOUBLE,MPI_MAX,_comm);
		if(timeArray[stopLevel]<minTime){
			minTime = timeArray[stopLevel];
			bestStopLevel = stopLevel;
		}
		int myRank;
		MPI_Comm_rank(_comm,&myRank);
		if(myRank == 0){
			std::cout << "stop level = " << stopLevel << " time = " << timeArray[stopLevel] << "\n";
		}

	}
	return bestStopLevel;
#else
	return 0;
#endif
}

void UniformPseudoParticleContainer::upwardPass(P2MCellProcessor* cp) {
	// P2M
	_leafContainer->traverseCells(*cp);
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL");

	int curCellsEdge=_globalNumCellsPerDim;
	double cellWid[3];
	if(_globalLevel >= _maxLevel - 1){
		communicateHalos();
	}
	for(int i=0; i <3; i++)	cellWid[i]=_cellLength[i];

	// when considering periodic boundary conditions, there is actually work up to level 1!
	for(int curLevel=_maxLevel-1; curLevel>=1; curLevel--){

		curCellsEdge /=2;
		for(int i=0; i <3; i++)	cellWid[i] *=2;
#if defined(ENABLE_MPI) && WIGNER==0
		if(curLevel >= _globalLevel){ //local M2M
			Vector3<int> curCellsEdgeLocal;
			for(int j = 0; j < 3; j++){
				curCellsEdgeLocal[j] = (int) (curCellsEdge/_numProcessorsPerDim[j])+4;
			}
			const Vector3<int> offset = (_globalLevel == curLevel)? _processorPositionGlobalLevel: Vector3<int>(2);
			CombineMpCell_Local(cellWid, curCellsEdgeLocal , curLevel, offset);
		}
		else{ //global M2M
			CombineMpCell_Global(cellWid, curCellsEdge, curLevel);
			if(_avoidAllReduce && curLevel >= _stopLevel){ //perform neighbourhood allreduce instead of global allreduce
				int mpCells = pow(2,curLevel);
				int stride = pow(2,_globalLevel - curLevel);
				int myRank;
				int coords[3];
				int coordsLevel[3];
				MPI_Comm_rank(_comm,&myRank);
				MPI_Cart_coords(_comm, myRank, 3, coords);
				for(int d = 0; d < 3; d++){
					coordsLevel[d] = ((coords[d] * _numCellsOnGlobalLevel[d]) / (stride));
				}
				int cellIndex = ((coordsLevel[2]) * mpCells + coordsLevel[1]) * mpCells + coordsLevel[0];
				if(!_fuseGlobalCommunication){
					auto buffer = std::vector<double>(2*_expansionSize);
					MpCell & currentCell = _mpCellGlobalTop[curLevel][cellIndex];
					int index = 0;
					currentCell.multipole.writeValuesToMPIBuffer(buffer,index);
					MPI_Allreduce(MPI_IN_PLACE,buffer.data(),2 * _expansionSize,MPI_DOUBLE,MPI_SUM,_neighbourhoodComms[curLevel]);
					index = 0;
					currentCell.multipole.readValuesFromMPIBuffer(buffer,index);
				}
				else{ //get values of 8 values from previous level that contributed to the parent level in addition to parent value (9 values each 2 expansions)
					int coordsFlooredPreviousLevel[3];
					for(int d = 0; d < 3; d++){
						coordsFlooredPreviousLevel[d] = (((coords[d] * _numCellsOnGlobalLevel[d]) / (stride/2)) / 2) * 2;
					}
					auto buffer = std::vector<double>(18*_expansionSize);
					MpCell & currentCell = _mpCellGlobalTop[curLevel][cellIndex];
					int index = 0;
					currentCell.multipole.writeValuesToMPIBuffer(buffer,index);
					mpCells *= 2; //mpCells for previous level
					for(int z = 0; z < 2; z++){
						for(int y = 0; y < 2; y++){
							for(int x = 0; x < 2; x++){
								int cellIndexTemp = ((coordsFlooredPreviousLevel[2] + z) * mpCells + coordsFlooredPreviousLevel[1] + y) * mpCells + coordsFlooredPreviousLevel[0] + x;
								MpCell & currentCellPrevious = _mpCellGlobalTop[curLevel + 1][cellIndexTemp];
								currentCellPrevious.multipole.writeValuesToMPIBuffer(buffer,index);

							}
						}
					}
					MPI_Allreduce(MPI_IN_PLACE,buffer.data(),18 * _expansionSize,MPI_DOUBLE,MPI_SUM,_neighbourhoodComms[curLevel]);
					index = 0;
					MpCell & currentCell2 = _mpCellGlobalTop[curLevel][cellIndex];
					currentCell2.multipole.readValuesFromMPIBuffer(buffer,index);
					for(int z = 0; z < 2; z++){
						for(int y = 0; y < 2; y++){
							for(int x = 0; x < 2; x++){
								int cellIndexTemp = ((coordsFlooredPreviousLevel[2] + z) * mpCells + coordsFlooredPreviousLevel[1] + y) * mpCells + coordsFlooredPreviousLevel[0] + x;
								MpCell & currentCellPrevious2 = _mpCellGlobalTop[curLevel + 1][cellIndexTemp];
								currentCellPrevious2.multipole.readValuesFromMPIBuffer(buffer,index);

							}
						}
					}

				}
			}
		}
		if(curLevel == _globalLevel + 1){
			communicateHalos();
		}
		if(curLevel == _globalLevel){
			global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL");
			global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_LOKAL");

		}
#else
		CombineMpCell_Global(cellWid, curCellsEdge, curLevel);
#endif

	}
#ifdef ENABLE_MPI
	if(_avoidAllReduce and _fuseGlobalCommunication and _stopLevel <= _globalLevel){ //get remaining 7 cells at stop level
		int mpCells = pow(2, _stopLevel);
		int stride = pow(2,_globalLevel - _stopLevel);
		int myRank;
		int coords[3];
		int coordsFlooredLevel[3];
		MPI_Comm_rank(_comm,&myRank);
		MPI_Cart_coords(_comm, myRank, 3, coords);
		for(int d = 0; d < 3; d++){
			coordsFlooredLevel[d] = (((coords[d] * _numCellsOnGlobalLevel[d]) / (stride)) / 2) * 2;
		}

		auto buffer = std::vector<double>(16*_expansionSize);
		int index = 0;
		for(int z = 0; z < 2; z++){
			for(int y = 0; y < 2; y++){
				for(int x = 0; x < 2; x++){
					int cellIndexTemp = ((coordsFlooredLevel[2] + z) * mpCells + coordsFlooredLevel[1] + y) * mpCells + coordsFlooredLevel[0] + x;
					MpCell & currentCell = _mpCellGlobalTop[_stopLevel][cellIndexTemp];
					currentCell.multipole.writeValuesToMPIBuffer(buffer,index);

				}
			}
		}
		MPI_Allreduce(MPI_IN_PLACE,buffer.data(),16 * _expansionSize,MPI_DOUBLE,MPI_SUM,_neighbourhoodComms[_stopLevel - 1]);
		index = 0;

		for(int z = 0; z < 2; z++){
			for(int y = 0; y < 2; y++){
				for(int x = 0; x < 2; x++){
					int cellIndexTemp = ((coordsFlooredLevel[2] + z) * mpCells + coordsFlooredLevel[1] + y) * mpCells + coordsFlooredLevel[0] + x;
					MpCell & currentCell = _mpCellGlobalTop[_stopLevel][cellIndexTemp];
					currentCell.multipole.readValuesFromMPIBuffer(buffer,index);

				}
			}
		}
	}
	if(_avoidAllReduce){
		if(_globalLevel >= 1){
			communicateOwnGlobalValue(_stopLevel);
			AllReduceMultipoleMomentsLevelToTop(pow(8,_stopLevel - 1), _stopLevel - 1);
		}
	}
	else{
		AllReduceMultipoleMomentsLevelToTop(_globalLevelNumCells,_globalLevel);
	}
	if(_globalLevel != 0){
		global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_LOKAL");
	}
	else{
		global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL");
	}
#else
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL");
#endif
//	//trigger communication by testing with MPI test
//	initBusyWaiting();
//	busyWaiting();

}

void UniformPseudoParticleContainer::horizontalPass(
		VectorizedChargeP2PCellProcessor* cp) {

	// P2P
	_leafContainer->traverseCellPairs(*cp);

	// M2L
	int curCellsEdge=1;
	double cellWid[3];

	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

#if defined(ENABLE_MPI)
	//iterate through local areas of local tree without halos
	for(int curLevel= 1; curLevel<=_maxLevel; curLevel++){

		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /=2;
		}
		if(curLevel > _globalLevel){
			Vector3<int> curCellsEdgeLocal;
			for(int j = 0; j < 3; j++){
				curCellsEdgeLocal[j] = (int) (curCellsEdge/_numProcessorsPerDim[j])+4;
			}
#ifdef FMM_FFT
			GatherWellSepLo_FFT_Local(cellWid, curCellsEdgeLocal, curLevel, 0);
#else
			GatherWellSepLo_Local(cellWid, curCellsEdgeLocal, curLevel, 0);
#endif
		}
	}
#endif

	int finishedFlag = 0;
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING");
	//initialize busy waiting variables
	initBusyWaiting();
	//start busy waiting
	while(finishedFlag != -1){
		finishedFlag = busyWaiting(); //returns -1 if finished
#if defined(ENABLE_MPI)
		if(finishedFlag == 1){ //communication of local tree halos finished
			global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING");
			communicateHalosOverlapSetHalos();
			//start receiving for second halo exchange in NT method
			if(_doNTLocal){
				_multipoleRecBufferOverlap->communicate(true);
			}
			curCellsEdge=1;
			for(int i=0; i < 3; i++) cellWid[i] = _domain->getGlobalLength(i);

			for(int curLevel=1; curLevel<=_maxLevel; curLevel++){ //do local tree halo M2M

				curCellsEdge *=2;
				for(int i=0; i < 3; i++){
					cellWid[i] /=2;
				}

				if(curLevel > _globalLevel){
					 Vector3<int> curCellsEdgeLocal;
					for(int j = 0; j < 3; j++){
						curCellsEdgeLocal[j] = (int) (curCellsEdge/_numProcessorsPerDim[j])+4;
					}
#ifdef FMM_FFT
					GatherWellSepLo_FFT_Local(cellWid, curCellsEdgeLocal, curLevel, 1);
#else
					GatherWellSepLo_Local(cellWid, curCellsEdgeLocal, curLevel, 1);
#endif
				}
			}
			global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING");
		}
#endif

		if(finishedFlag == 2){ //communication of global tree finished
			global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING");
#ifdef ENABLE_MPI
			if(_avoidAllReduce){ // set halo values
				if(_globalLevel >= 2){
					AllReduceMultipoleMomentsSetValues(pow(8,_stopLevel - 1), _stopLevel - 1);
				}
				communicateHaloGlobalValues(_stopLevel);
			}
			else{
				AllReduceMultipoleMomentsSetValues(pow(8,_globalLevel), _globalLevel);
			}
			if(_doNTGlobal and _avoidAllReduce and not (_globalLevel == 1 and _fuseGlobalCommunication)){
				_multipoleRecBufferOverlapGlobal->communicateGlobalLevels(_globalLevel,_stopLevel,true);
			}
#endif
			curCellsEdge=1;
			for(int i=0; i < 3; i++) cellWid[i] = _domain->getGlobalLength(i);

			for(int curLevel=1; curLevel<=_maxLevel; curLevel++){ //global M2M

				curCellsEdge *=2;
				for(int i=0; i < 3; i++){
					cellWid[i] /=2;
				}

#if defined(ENABLE_MPI)
				if(curLevel <= _globalLevel){
	#ifdef FMM_FFT
					GatherWellSepLo_FFT_Global(cellWid, curCellsEdge, curLevel);
	#else
					GatherWellSepLo_Global(cellWid, curCellsEdge, curLevel);
	#endif
				}
#else
	#ifdef FMM_FFT
					GatherWellSepLo_FFT_Global(cellWid, curCellsEdge, curLevel);
	#else
					GatherWellSepLo_Global(cellWid, curCellsEdge, curLevel);
	#endif
					finishedFlag = -1;
#endif
#ifdef ENABLE_MPI
				if(_doNTGlobal && _avoidAllReduce && _fuseGlobalCommunication && curLevel >= _stopLevel && curLevel <= _globalLevel){ //perform local allreduce for backwards communication
					int mpCells = pow(2,curLevel);
					int stride = pow(2,_globalLevel - curLevel);
					int myRank;
					int coords[3];
					int coordsLevel[3];
					MPI_Comm_rank(_comm,&myRank);
					MPI_Cart_coords(_comm, myRank, 3, coords);
					//get values of 216 halo values or 8 if level == 1
					int coordsFlooredLevel[3];
					for(int d = 0; d < 3; d++){
						coordsFlooredLevel[d] = (((coords[d] * _numCellsOnGlobalLevel[d]) / stride) / 2) * 2;
					}
					//Fixme? special case for curlevel == 2? -> 64 cells only? unnecessary communication in this case?
					//possible optimization only add up values that were modified in NT method
					int numCells = (curLevel == 1)? 8 : 216;
					auto buffer = std::vector<double>(numCells * 2 * _expansionSize);
					int index = 0;
					int start, end;
					if(curLevel == 1){
						start = 0;
						end = 2;
					}
					else{
						start = -2;
						end = 4;
					}
					//add up all the halo contributions of the 8 processors
					for(int z = start; z < end; z++){
						for(int y = start; y < end; y++){
							for(int x = start; x < end; x++){
								int zIndex = (coordsFlooredLevel[2] + z + mpCells) % mpCells;
								int yIndex = (coordsFlooredLevel[1] + y + mpCells) % mpCells;
								int xIndex = (coordsFlooredLevel[0] + x + mpCells) % mpCells;

								int cellIndexTemp = (zIndex * mpCells + yIndex) * mpCells + xIndex;
								MpCell & currentCell = _mpCellGlobalTop[curLevel][cellIndexTemp];
								currentCell.local.writeValuesToMPIBuffer(buffer,index);

							}
						}
					}
					MPI_Allreduce(MPI_IN_PLACE,buffer.data(),numCells * 2 * _expansionSize,MPI_DOUBLE,MPI_SUM,_neighbourhoodComms[curLevel - 1]);
					index = 0;

					for(int z = start; z < end; z++){
						for(int y = start; y < end; y++){
							for(int x = start; x < end; x++){
								int zIndex = (coordsFlooredLevel[2] + z + mpCells) % mpCells;
								int yIndex = (coordsFlooredLevel[1] + y + mpCells) % mpCells;
								int xIndex = (coordsFlooredLevel[0] + x + mpCells) % mpCells;

								int cellIndexTemp = (zIndex * mpCells + yIndex) * mpCells + xIndex;
								MpCell & currentCell = _mpCellGlobalTop[curLevel][cellIndexTemp];
								currentCell.local.readValuesFromMPIBuffer(buffer,index);

							}
						}
					}

				}
#endif
			}
			global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING");
		}
#if defined(ENABLE_MPI)
		if(finishedFlag == 5){ //local back communication receive finished -> set values
			communicateHalosOverlapPostProcessingSetHalos();
		}
		if(finishedFlag == 4){ //local halos processed and sending and receiving finished -> start back communication
			communicateHalosOverlapPostProcessingStart();
		}

		if(finishedFlag == 7){ //global back communication receive finished -> set values
			if(_avoidAllReduce){
				if(_globalLevel >= 1){
					communicateOwnGlobalValue(_stopLevel,true);
				}
			}
		}
		if(finishedFlag == 6){ //global halos processed and sending and receiving finished -> start back communication
			if(_avoidAllReduce){
				if(_globalLevel >= 1){
					communicateHaloGlobalValues(_stopLevel,true);
				}
			}
		}
#endif
	}
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING");
	//in case of neutral territory version exchange halo values
}

int UniformPseudoParticleContainer::busyWaiting(){
#ifdef ENABLE_MPI
	//if program has dead lock comment this line in to see which flag causes dead lock
//	std::cout << _allReduceProcessed << _halosProcessed << _sendLocalProcessed << _sendGlobalProcessed << _backCommunicationLocalProcessed << _backCommunicationGlobalProcessed << _globalHalosProcessed <<"\n";
	if(_allReduceProcessed and _halosProcessed and _sendLocalProcessed and _sendGlobalProcessed and _backCommunicationLocalProcessed and _backCommunicationGlobalProcessed and _globalHalosProcessed){
		return -1;
	}
	if(_halosProcessed and _sendLocalProcessed and !_backCommunicationLocalStarted){
		if(_doNTLocal){
			_backCommunicationLocalStarted = 1;
			return 4;
		}
	}
	if(!_backCommunicationGlobalProcessed and _globalHalosProcessed and _allReduceProcessed and _sendGlobalProcessed and !_backCommunicationGlobalStarted){
		if(_doNTGlobal){
			_backCommunicationGlobalStarted = 1;
			return 6;
		}
	}
	if((!_backCommunicationLocalProcessed) and _backCommunicationLocalStarted){

		if(_multipoleBufferOverlap->testIfFinished() and _multipoleRecBufferOverlap->testIfFinished()){
			_backCommunicationLocalProcessed = 1;
			return 5;
		}
	}
	if((!_backCommunicationGlobalProcessed) and _backCommunicationGlobalStarted){

		if(_multipoleBufferOverlapGlobal->testIfFinished() and _multipoleRecBufferOverlapGlobal->testIfFinished()){
			_backCommunicationGlobalProcessed = 1;
			return 7;
		}
	}
	if(!_sendLocalProcessed){
		if(_multipoleBufferOverlap->testIfFinished()){
			_sendLocalProcessed = 1;

			return 3;
		}
	}
	if(!_sendGlobalProcessed){
		if(_multipoleBufferOverlapGlobal->testIfFinished()){
			_sendGlobalProcessed = 1;

			return 8;
		}
	}
	if(!_halosProcessed){
		if(_multipoleRecBufferOverlap->testIfFinished()){
			_halosProcessed = 1;

			return 1;
		}
	}
	if(!_globalHalosProcessed){
		if(_multipoleRecBufferOverlapGlobal->testIfFinished()){
			_globalHalosProcessed = 1;
			if(_allReduceProcessed){
				return 2;
			}
		}
	}
	if(!_allReduceProcessed){
		int flag;
		MPI_Status status;
		MPI_Test(&_allReduceRequest, &flag, &status);
		if(flag){
			_allReduceProcessed = 1;
			if(_globalHalosProcessed){
				return 2;
			}
		}
	}
	return 0;
#else
	return 2;
#endif
}

void UniformPseudoParticleContainer::downwardPass(L2PCellProcessor* cp) {
	// L2L
	int curCellsEdge=1;
	double cellWid[3];
#ifdef ENABLE_MPI
	//start receiving for next iteration; important for ready send
	_multipoleRecBufferOverlap->communicate(false);
	if(_avoidAllReduce){
		//start receiving for next iteration
		if(_globalLevel >= 1 and not (_globalLevel == 1 and _fuseGlobalCommunication)){
			_multipoleRecBufferOverlapGlobal->communicateGlobalLevels(_globalLevel,_stopLevel);
		}
	}
#endif
	for(int i=0; i <3; i++)
		cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<_maxLevel; curLevel++){
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /= 2;
		}

#if defined(ENABLE_MPI) && WIGNER==0
		if(curLevel >= _globalLevel){ //do local M2M
			Vector3<int> curCellsEdgeLocal;
			for(int j = 0; j < 3; j++){
				curCellsEdgeLocal[j] = (int) (curCellsEdge/_numProcessorsPerDim[j])+4;
			}
			const Vector3<int> offset = (_globalLevel == curLevel)? _processorPositionGlobalLevel: Vector3<int>(2);
			PropagateCellLo_Local(cellWid, curCellsEdgeLocal, curLevel,offset);
		}
		else{ //do global M2M
			PropagateCellLo_Global(cellWid, curCellsEdge, curLevel);
		}
#else
		PropagateCellLo_Global(cellWid, curCellsEdge, curLevel);
#endif
	}

	// L2P
	_leafContainer->traverseCells(*cp);

	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_FMM_COMPLETE");
}

void UniformPseudoParticleContainer::CombineMpCell_Global(double */*cellWid*/, int mpCells, int curLevel){
	int iDir,
		m1       = 0,
		m1x,
		m1y,
		m1z,
		m2       = 0;
	int m2v[3]   = {0, 0, 0};
	int mpCellsN = 2 * mpCells;

	for (int mloop    = 0; mloop < mpCells * mpCells * mpCells; mloop++) {
		m1x = mloop % mpCells;
		m1y = (mloop / mpCells) % mpCells;
		m1z = mloop / (mpCells * mpCells);
		m1  = (m1z * mpCells + m1y) * mpCells + m1x;

		for (iDir = 0; iDir < 8; iDir++) { //iterate over 8 children of m1

			m2v[0] = 2 * m1x;
			m2v[1] = 2 * m1y;
			m2v[2] = 2 * m1z;

			if (IsOdd(iDir)) m2v[0]     = m2v[0] + 1;
			if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
			if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;


			m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];

			if (_mpCellGlobalTop[curLevel + 1][m2].occ == 0) continue;

			_mpCellGlobalTop[curLevel][m1].occ += _mpCellGlobalTop[curLevel + 1][m2].occ;

			_mpCellGlobalTop[curLevel][m1].multipole.addMultipoleParticle(_mpCellGlobalTop[curLevel + 1][m2].multipole);
		} // iDir closed
	} // mloop closed
}

void UniformPseudoParticleContainer::CombineMpCell_Local(double* /*cellWid*/, Vector3<int> localMpCells, int curLevel, Vector3<int> offset){
	int iDir,
		m1     = 0,
		m1x,
		m1y,
		m1z,
		m2     = 0;
	int m2v[3] = {0, 0, 0};
	//take care of halo cells
	Vector3<int> localMpCellsN;
	for(int i = 0; i < 3; ++i){
		localMpCellsN[i] = 2 * (localMpCells[i] - 4) + 4;
	}
	Vector3<int> localMpCellsRow;
	std::vector<std::vector<MpCell> > * mpCellCurLevel;
	int curLevelp1 = curLevel -_globalLevel;

	if(curLevel == _globalLevel){
		for(int i = 0; i < 3; ++i){
			localMpCellsRow[i] = (localMpCells[i]-4) * _numProcessorsPerDim[i];
		}
		mpCellCurLevel = &_mpCellGlobalTop;

	}
	else{
		localMpCellsRow = localMpCells;
		mpCellCurLevel = &_mpCellLocal;
		//adjust level to local tree
		curLevel = curLevel - _globalLevel - 1;
	}
	//current level plus 1
	int numInnerCells[3];
	for(int d = 0; d < 3; ++d){
		numInnerCells[d] = localMpCells[d] - 4;
	}

	for (int mloop = 0 ; mloop < numInnerCells[0] * numInnerCells[1] * numInnerCells[2]; mloop++){
		m1x = mloop % numInnerCells[0];
		m1y = (mloop / numInnerCells[0]) % numInnerCells[1];
		m1z = mloop / (numInnerCells[0] * numInnerCells[1]);
		m1=((m1z+offset[2])*localMpCellsRow[1] + m1y+offset[1])*localMpCellsRow[0] + m1x+offset[0];

		for(iDir=0; iDir<8; ++iDir){ //iterate over children

			m2v[0]=2*m1x+2;
			m2v[1]=2*m1y+2;
			m2v[2]=2*m1z+2;

			if(IsOdd(iDir  )) m2v[0]=m2v[0]+1;
			if(IsOdd(iDir/2)) m2v[1]=m2v[1]+1;
			if(IsOdd(iDir/4)) m2v[2]=m2v[2]+1;


			m2=(m2v[2]*localMpCellsN[1] + m2v[1])*localMpCellsN[0] + m2v[0];

			if(_mpCellLocal[curLevelp1][m2].occ==0) continue;

			(*mpCellCurLevel)[curLevel][m1].occ +=_mpCellLocal[curLevelp1][m2].occ;
			(*mpCellCurLevel)[curLevel][m1].multipole.addMultipoleParticle(_mpCellLocal[curLevelp1][m2].multipole);

		} // iDir closed
		if(curLevel == _globalLevel){
			(*mpCellCurLevel)[curLevel][m1].occ++; //ensure that complete subtree is never considered to be empty
		}
	} // mloop closed
}

#define HiLim(t) ToEven(m1v[t])+ 2*_wellSep+1
#define LoLim(t) ToEven(m1v[t])- 2*_wellSep

void UniformPseudoParticleContainer::GatherWellSepLo_Global(double *cellWid, int mpCells, int curLevel){
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL");
	int m1v[3];
	int m2v[3];
	int m1,
		m2,
		m2x,
		m2y,
		m2z;
	// int m1x, m1y, m1z;
	int m22x,
		m22y,
		m22z; // for periodic image
	int _row_length;
	_row_length = mpCells * mpCells * mpCells;
	DomainDecompBase &domainDecomp = global_simulation->domainDecomposition();
	Vector3<double> periodicShift;
	for (int        m1Loop         = 0; m1Loop < _row_length; m1Loop++) {

		m1v[0] = m1Loop % mpCells;
		m1v[1] = (m1Loop / mpCells) % mpCells;
		m1v[2] = (m1Loop / (mpCells * mpCells)) % mpCells;
		if (_mpCellGlobalTop[curLevel][m1Loop].occ == 0)
			continue;
		int offsetEnd, offsetStart;
		if (!_doNTGlobal or curLevel < _stopLevel) { // no NT in this case
			offsetEnd   = 0;
			offsetStart = 0;
		} else { //iterate over tower if NT case
			offsetStart = -2;
			offsetEnd   = 3;
		}
		int      m1v1_local = m1v[1];
		for (int yOffset    = offsetStart; yOffset <= offsetEnd; yOffset++) { //iterate over tower

			if (offsetEnd != 0 or offsetStart != 0) {
				m1v[1] = (m1v1_local & ~1) + yOffset;
				int m11y = ((m1v1_local & ~1) + yOffset + mpCells) % mpCells;
				m1 = (m1v[2] * mpCells + m11y) * mpCells + m1v[0];
			} else {
				m1 = (m1v[2] * mpCells + m1v[1]) * mpCells + m1v[0];
			}
			for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
				if (_periodicBC == false and (m2z < 0 or m2z >= mpCells)) {
					continue;
				}
				// to get periodic image
				m22z = (mpCells + m2z) % mpCells;
				periodicShift[2]                     = 0.0;
				if (m2z < 0) periodicShift[2]        = -mpCells * cellWid[2];
				if (m2z >= mpCells) periodicShift[2] = mpCells * cellWid[2];

				m2v[2]   = m2z;
				int m2yStart, m2yEnd;
				if (_doNTGlobal && curLevel >= _stopLevel) { //guarantees that m2y stays in y interval of plate
					m2yStart = m2yEnd = m1v1_local;
				} else {
					m2yStart = LoLim(1);
					m2yEnd   = HiLim(1);
				}
				for (m2y = m2yStart; m2y <= m2yEnd; m2y++) {
					if (_periodicBC == false and (m2y < 0 or m2y >= mpCells)) {
						continue;
					}
					// to get periodic image
					m22y = (mpCells + m2y) % mpCells;
					periodicShift[1] = 0.0;
					if (_doNTGlobal && curLevel >= _stopLevel) { //guarantees that m2y stays in y interval of plate
						if (m1v[1] < 0) periodicShift[1]        = mpCells * cellWid[1];
						if (m1v[1] >= mpCells) periodicShift[1] = -mpCells * cellWid[1];
					} else {
						if (m2y < 0) periodicShift[1]        = -mpCells * cellWid[1];
						if (m2y >= mpCells) periodicShift[1] = mpCells * cellWid[1];
					}
					m2v[1]           = m2y;
					for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
						if (_periodicBC == false and (m2x < 0 or m2x >= mpCells)) {
							continue;
						}
						// to get periodic image
						m22x = (mpCells + m2x) % mpCells;
						periodicShift[0]                     = 0.0;
						if (m2x < 0) periodicShift[0]        = -mpCells * cellWid[0];
						if (m2x >= mpCells) periodicShift[0] = mpCells * cellWid[0];
						//
						m2v[0]                               = m2x;
						m2 = (m22z * mpCells + m22y) * mpCells + m22x;
						if (filterM2Global(curLevel, m2v, m1v, m2x, m2y, m2z, m2, yOffset)) {
							continue;
						}
						_mpCellGlobalTop[curLevel][m1].local.addMultipoleParticle(
								_mpCellGlobalTop[curLevel][m2].multipole, periodicShift);
						if (_doNTGlobal && curLevel >= _stopLevel) { //for NT do both directions
							_mpCellGlobalTop[curLevel][m2].local.addMultipoleParticle(
									_mpCellGlobalTop[curLevel][m1].multipole, -1 * periodicShift);
						}
					} // m2x closed
				} // m2y closed
			} // m2z closed
		} // tower closed
	} //m1 closed
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL");
} // GatherWellSepLo closed

bool UniformPseudoParticleContainer::filterM2Global(int curLevel, int *m2v, int *m1v, int m2x, int m2y, int m2z, int m2, int yOffset){
	//check if well separated
	if (abs(m2v[0] - m1v[0]) <= _wellSep &&
		abs(m2v[1] - m1v[1]) <= _wellSep &&
		abs(m2v[2] - m1v[2]) <= _wellSep){
		return true;
	}
	//only accept if m2 is in plate and avoid double calculation of forces
	if(_doNTGlobal && curLevel >= _stopLevel and (not(m2x >= LoLim(0) + 2 and not(m2x <= HiLim(0) - 2 and m2z < LoLim(2) + 2)) //plate
	   or (yOffset < 0 and m2x <= HiLim(0) - 2 and m2z <= HiLim(2) - 2))){ // no interaction from local area to lower tower
		return true;
	}
	//avoid empty cells; not used in MPI version as cells might be non-empty in global tree even with occ=0
#ifndef ENABLE_MPI
	if (_mpCellGlobalTop[curLevel][m2].occ == 0)
		return true;
#endif
	return false;
}
bool UniformPseudoParticleContainer::filterM1Local(bool doHalos, int m1, int m1x, int m1y, int m1z, Vector3<int> localMpCells, int curLevel){
	//accept only border cells as all other inner computations have been completed in previous call with doHalos=false
	if(doHalos and (m1x >= 4 and m1y >= 4 and m1z >= 4 and m1x < localMpCells[0] - 4 and m1y < localMpCells[1] - 4 and m1z < localMpCells[2] -4)
			//or (m1z < 2 or m1z >= localMpCells[2] - 2 or m1y < 2 or m1y >= localMpCells[1] - 2 or m1x < 2 or m1x >= localMpCells[0] - 2)) //cannot be halo cell
			){ // either inner or halo cell
		return true;
	}
	//accept only if m1 is in tower when using NT method
	if(doHalos and _doNTLocal and (m1x < 2 or m1z < 2 or m1x >= localMpCells[0] - 2 or m1z >= localMpCells[2] - 2)){ //m1 should be in local subregion or in adjacent area region in y direction (tower)
		return true;
	}
	//do not accept empty cells for m1
	if (_mpCellLocal[curLevel][m1].occ == 0 ){
		return true;
	}
	return false;
}

bool UniformPseudoParticleContainer::filterM2Local(bool doHalos, int m1, int m1x, int m1y, int m1z, int m2, int m2x, int m2y, int m2z, Vector3<int> localMpCells, int curLevel, bool inHaloz, bool inHaloy, bool inHalox){
	if(not(_doNTLocal)){
		if(doHalos and not(inHalox or inHaloy or inHaloz)){ //during halo calculation one of the cells needs to be in the halo region
			return true;
		}
	}
	if(doHalos and _doNTLocal and (m2x < 2 or (m2x < localMpCells[0] - 2 and m2z < 2) or //m2 must be in second NT region
			(m1y < 2 and not(inHalox) and not(inHaloz)) or  //local subregion may not interact with lower y half
			(m1y >= 2 and m1y < localMpCells[1] - 2 and not(inHalox) and not(inHaloz)))){  //not both may be in local subregion (those interactions were handled in previous call with doHalos=0)
//		inHalox = 0; // must be done in calling function
		return true;
	}
	//check if well separated
	if (abs(m2x - m1x) <= _wellSep &&
		abs(m2y - m1y) <= _wellSep &&
		abs(m2z - m1z) <= _wellSep){
		return true;
	}
	//check if cell occupied
	if (_mpCellLocal[curLevel][m2].occ == 0){
		return true;
	}
	return false;
}

void UniformPseudoParticleContainer::GatherWellSepLo_Local(double* /*cellWid*/, Vector3<int> localMpCells, int curLevel, int doHalos){
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_LOKAL");
	if (doHalos) {
		global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_HALO_GATHER");
	}
	int m1x, m1y, m1z;
	int m1v[3];
	//adjust for local level
	curLevel = curLevel - _globalLevel - 1;
	int             m1, m2, m2x, m2y, m2z;
	int             _row_length;
	Vector3<double> periodicShift(0.0);
	int             offset;
	if (_doNTLocal) {
		offset = 0;
	} else {
		offset = 2;
	}
	int      zStart = offset;
	int      xStart = offset;
	int      yStart = offset;
	int      xEnd   = localMpCells[0] - offset - xStart;
	int      yEnd   = localMpCells[1] - offset - yStart;
	int      zEnd   = localMpCells[2] - offset - zStart;
	for (int mloop  = 0; mloop < xEnd * yEnd * zEnd; mloop++) {
		m1x = mloop % xEnd + xStart;
		m1y = (mloop / xEnd) % yEnd + yStart;
		m1z = mloop / (xEnd * yEnd) + zStart;
		m1 = ((m1z) * localMpCells[1] + m1y) * localMpCells[0] + m1x;
		if (filterM1Local(doHalos, m1, m1x, m1y, m1z, localMpCells, curLevel)) { //check if m1 should be skipped
			continue;
		}
		m1v[0] = m1x;
		m1v[1] = m1y;
		m1v[2] = m1z;
		bool inHaloz, inHaloy, inHalox; //shows if current cell is in halo area in respective coordinate axis
		inHaloz  = inHaloy = inHalox = 0;
		for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
			if ((m2z < 2 or m2z >= localMpCells[2] - 2)) { //halo cell
				if (!doHalos) {
					continue;
				}
				if (doHalos) {
					inHaloz = 1;
				}
			}
			for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
				if ((m2y < 2 or m2y >= localMpCells[1] - 2)) { //halo cell
					if (!doHalos) {
						continue;
					}
					if (doHalos) {
						if (_doNTLocal) {//m2 not in y halo allowed in NT (m2 in plate and not in tower)
							continue;
						} else {
							inHaloy = 1;
						}
					}
				}
				for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
					if ((m2x < 2 or m2x >= localMpCells[0] - 2)) { //halo cell
						if (!doHalos) {
							continue;
						}
						if (doHalos) {
							inHalox = 1;
						}
					}
					m2      = (m2z * localMpCells[1] + m2y) * localMpCells[0] + m2x;
					if (filterM2Local(doHalos, m1, m1x, m1y, m1z, m2, m2x, m2y, m2z, localMpCells, curLevel, inHaloz,
									  inHaloy, inHalox)) { //check if this m2 value needs to be skipped
						inHalox = 0;
						continue;
					}
					inHalox = 0;
					_mpCellLocal[curLevel][m1].local.addMultipoleParticle(
							_mpCellLocal[curLevel][m2].multipole, periodicShift);
					if (_doNTLocal and doHalos) {
						_mpCellLocal[curLevel][m2].local.addMultipoleParticle(
								_mpCellLocal[curLevel][m1].multipole, -1 * periodicShift);
					}
				} // m2x closed
				inHaloy  = 0;
			} // m2y closed
			inHaloz  = 0;
		} // m2z closed

	} //mloop closed
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_LOKAL");
	if (doHalos) {
		global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_HALO_GATHER");
	}
} // GatherWellSepLo closed

#ifdef FMM_FFT
void UniformPseudoParticleContainer::GatherWellSepLo_FFT_Global(double *cellWid, int mpCells, int curLevel) {
	if (FFTSettings::USE_VECTORIZATION) {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Global_template<true, true, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_Global_template<true, true, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Global_template<true, true, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_Global_template<true, true, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Global_template<true, false, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_Global_template<true, false, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Global_template<true, false, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_Global_template<true, false, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		}
	} else {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Global_template<false, true, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_Global_template<false, true, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Global_template<false, true, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_Global_template<false, true, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Global_template<false, false, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_Global_template<false, false, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Global_template<false, false, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_Global_template<false, false, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		}
	}
}

void UniformPseudoParticleContainer::GatherWellSepLo_FFT_Local(double *cellWid, Vector3<int> mpCells, int curLevel, int doHalos) {
	if (FFTSettings::USE_VECTORIZATION) {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Local_template<true, true, true, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_Local_template<true, true, true, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Local_template<true, true, false, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_Local_template<true, true, false, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Local_template<true, false, true, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_Local_template<true, false, true, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Local_template<true, false, false, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_Local_template<true, false, false, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			}
		}
	} else {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Local_template<false, true, true, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_Local_template<false, true, true, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Local_template<false, true, false, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_Local_template<false, true, false, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Local_template<false, false, true, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_Local_template<false, false, true, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_Local_template<false, false, false, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_Local_template<false, false, false, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			}
		}
	}
}

template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
void UniformPseudoParticleContainer::GatherWellSepLo_FFT_Global_template(
		double *cellWid, int mpCells, int curLevel) {
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL");
	int m1v[3]; //cell1 coordinates
	int m2v[3]; //cell2 coordinates
	int m1, m2, m2x, m2y, m2z; //cell coordinates
	int m22x, m22y, m22z; // for periodic image
	int loop_min, loop_max;
	int row_length; //number of cells per row
	row_length = mpCells * mpCells * mpCells;
	loop_min   = 0;
	loop_max   = row_length;
	//FFT param
	double radius;
	double base_unit = 2.0 / sqrt(3);
	int M2L_order;
	FFTDataContainer* tf;
	//Initialize FFT
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_INIT");
	// FFT Initialize
	for (m1 = loop_min; m1 < loop_max; m1++) {
#ifdef ENABLE_MPI
		if(curLevel > 2 or (_doNTGlobal and curLevel > 1)){ //for low levels it is better to initialize all cells directly to avoid double initialization
#endif
		if (_mpCellGlobalTop[curLevel][m1].occ == 0)
			continue;
#ifdef ENABLE_MPI
		}
#endif
		radius = _mpCellGlobalTop[curLevel][m1].local.getRadius();
		FFTAccelerableExpansion &source =
										static_cast<bhfmm::SHMultipoleParticle &>(_mpCellGlobalTop[curLevel][m1].multipole).getExpansion();
		FFTAccelerableExpansion &target =
										static_cast<bhfmm::SHLocalParticle &>(_mpCellGlobalTop[curLevel][m1].local).getExpansion();
		_FFTAcceleration->FFT_initialize_Source(source, radius);
		_FFTAcceleration->FFT_initialize_Target(target);
#ifdef ENABLE_MPI
		if(curLevel > 2 or (_doNTGlobal and curLevel > 1)){ //initialize only cells that are accessed later during M2M calculation -> same traversal scheme
			//iterate over all accessible cells m2 (they are always occ=0 -> extra treatment)
			m1v[0] = m1 % mpCells;
			m1v[1] = (m1 / mpCells) % mpCells;
			m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
			if(_numCellsOnGlobalLevel[1] == 2 && m1v[1] % 2 == 1){ //don't iterate over tower twice
				continue;
			}
			for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
				// to get periodic image
				m22z = (mpCells + m2z) % mpCells;

				m2v[2] = m2z;
				for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
					// to get periodic image
					m22y = (mpCells + m2y) % mpCells;
					m2v[1] = m2y;
					for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
						// to get periodic image
						m22x = (mpCells + m2x) % mpCells;
						if(_doNTGlobal && curLevel >= _stopLevel and not(m2x == m1v[0] and m2z == m1v[2]) //tower
							and not(((_numCellsOnGlobalLevel[1] == 2 and m2y <= HiLim(1) - 2 and m2y >= LoLim(1) + 2) or m1v[1] == m2y)  and m2x >= LoLim(0) + 2 and not(m2x <= HiLim(0) - 2 and m2z < LoLim(2) + 2))){ //plate
							continue; //if not in tower or plate
						}
						m2v[0] = m2x;
						m2 = (m22z * mpCells + m22y) * mpCells + m22x;

						if(_doNTGlobal and curLevel == 2 and ( m2y < LoLim(1) + 2 or (m2x > HiLim(0) - 2 and m2z < LoLim(2) + 2))){ //avoid double initialization in case of level == 2
							continue;
						}
						//
						if(!_doNTGlobal  or curLevel < _stopLevel){ //check if well-separated
							if (abs(m2v[0] - m1v[0]) <= _wellSep
									&& abs(m2v[1] - m1v[1]) <= _wellSep
									&& abs(m2v[2] - m1v[2]) <= _wellSep)
								continue;
						}
						radius = _mpCellGlobalTop[curLevel][m2].local.getRadius();
						FFTAccelerableExpansion& source2 =
								static_cast<bhfmm::SHMultipoleParticle&>(_mpCellGlobalTop[curLevel][m2].multipole).getExpansion();
						FFTAccelerableExpansion& target2 =
								static_cast<bhfmm::SHLocalParticle&>(_mpCellGlobalTop[curLevel][m2].local).getExpansion();
						_FFTAcceleration->FFT_initialize_Source(source2, radius);
						_FFTAcceleration->FFT_initialize_Target(target2);
					}
				}
			}
		}
#endif
	}
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_INIT");
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL");
	//M2L in Fourier space
	for (int m1Loop = loop_min; m1Loop < loop_max; m1Loop++) {
		m1v[0] = m1Loop % mpCells;
		m1v[1] = (m1Loop / mpCells) % mpCells;
		m1v[2] = (m1Loop / (mpCells * mpCells)) % mpCells;
		if (_mpCellGlobalTop[curLevel][m1Loop].occ == 0)
			return;
		int offsetEnd, offsetStart;
		if(!_doNTGlobal or curLevel < _stopLevel){
			offsetEnd = 0;
			offsetStart = 0;
		}
		else{ //iterate over tower
			offsetStart = -2;
			offsetEnd = 3;
		}
		int m1v1_local = m1v[1];
		for(int yOffset = offsetStart ; yOffset <=offsetEnd; yOffset++){ //iterate over tower

			if(offsetEnd != 0 or offsetStart != 0){
				m1v[1] = (m1v1_local & ~1) + yOffset;
				int m11y = ((m1v1_local & ~1) + yOffset + mpCells) % mpCells;
				m1 = (m1v[2] * mpCells + m11y) * mpCells + m1v[0];
			}
			else{
				m1 = (m1v[2] * mpCells + m1v[1]) * mpCells + m1v[0];

			}
			for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
				// to get periodic image
				m22z = (mpCells + m2z) % mpCells;

				m2v[2] = m2z;
				int m2yStart, m2yEnd;
				if(_doNTGlobal && curLevel >= _stopLevel){ //guarantees that m2y stays in y interval of plate
					m2yStart = m2yEnd = m1v1_local;
				}
				else{
					m2yStart = LoLim(1);
					m2yEnd = HiLim(1);
				}
				for (m2y = m2yStart; m2y <= m2yEnd; m2y++) {
					// to get periodic image
					m22y = (mpCells + m2y) % mpCells;

					m2v[1] = m2y;
					for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
						// to get periodic image
						m22x = (mpCells + m2x) % mpCells;
						m2v[0] = m2x;
						m2 = (m22z * mpCells + m22y) * mpCells + m22x;
						if(filterM2Global(curLevel, m2v, m1v, m2x, m2y, m2z, m2,yOffset)){
							continue;
						}
#ifndef ENABLE_MPI
						if (UseM2L_2way) {
							if (m1 > m2)
								continue;
						}
#endif
						int doBothDirections = (_doNTGlobal && curLevel >= _stopLevel and not UseM2L_2way)? 1 : 0;
						for(int i = 0; i <= doBothDirections; i++){
							tf = _FFT_TM->getTransferFunction(m2v[0] - m1v[0],
															  m2v[1] - m1v[1], m2v[2] - m1v[2], base_unit, base_unit,
															  base_unit);
#ifndef _OPENMP
						// Can only be used in sequential mode because race conditions
						global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL");
						global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_CALCULATION");
#endif
							if (UseM2L_2way) {
								FFTAccelerableExpansion& source2 =
															   static_cast<bhfmm::SHMultipoleParticle&>(_mpCellGlobalTop[curLevel][m2].multipole).getExpansion();
								FFTAccelerableExpansion& target2 =
															   static_cast<bhfmm::SHLocalParticle&>(_mpCellGlobalTop[curLevel][m2].local).getExpansion();
								FFTAccelerableExpansion& source1 =
															   static_cast<bhfmm::SHMultipoleParticle&>(_mpCellGlobalTop[curLevel][m1].multipole).getExpansion();
								FFTAccelerableExpansion& target1 =
															   static_cast<bhfmm::SHLocalParticle&>(_mpCellGlobalTop[curLevel][m1].local).getExpansion();

								if (UseOrderReduction) {
									M2L_order = FFTOrderReduction::getM2LOrder(
											m2v[0] - m1v[0], m2v[1] - m1v[1],
											m2v[2] - m1v[2], _maxOrd);
									if (UseVectorization) {
										static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_2way_ORed_vec(
												source2, source1, target2, target1, tf,
												M2L_order);
									} else {
										static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_2way_ORed(
												source2, source1, target2, target1, tf,
												M2L_order);
									}
								} else {
									if (UseVectorization) {
										static_cast<FFTAccelerationAPI_2Way*>(_FFTAcceleration)->FFT_M2L_2way_vec(
												source2, source1, target2, target1, tf);
									} else {
										static_cast<FFTAccelerationAPI_2Way*>(_FFTAcceleration)->FFT_M2L_2way(
												source2, source1, target2, target1, tf);
									}
								}
							} else {
								FFTAccelerableExpansion& source =
															   static_cast<bhfmm::SHMultipoleParticle&>(_mpCellGlobalTop[curLevel][m2].multipole).getExpansion();
								FFTAccelerableExpansion& target =
															   static_cast<bhfmm::SHLocalParticle&>(_mpCellGlobalTop[curLevel][m1].local).getExpansion();

								if (UseOrderReduction) {
									M2L_order = FFTOrderReduction::getM2LOrder(
											m2v[0] - m1v[0], m2v[1] - m1v[1],
											m2v[2] - m1v[2], _maxOrd);
									if (UseVectorization) {
										static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_OrderReduction_vec(
												source, target, tf, M2L_order);
									} else {
										static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_OrderReduction(
												source, target, tf, M2L_order);
									}
								} else {
									if (UseVectorization) {
										_FFTAcceleration->FFT_M2L_vec(source, target,
																	  tf);
									} else {
										_FFTAcceleration->FFT_M2L(source, target, tf);
									}
								}
							}
#ifndef _OPENMP
						// Can only be used in sequential mode because race conditions
						global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_CALCULATION");
						global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL");
#endif
							if (!UseTFMemoization) {
								delete tf; //free useless memory
							}
							if(_doNTGlobal && curLevel >= _stopLevel and not UseM2L_2way){
								//exchange m1 and m2 for second direction
								int temp = m1;
								m1 = m2;
								m2 = temp;
								int tempAr[3];
								tempAr[0] = m1v[0];
								tempAr[1] = m1v[1];
								tempAr[2] = m1v[2];
								m1v[0] = m2v[0];
								m1v[1] = m2v[1];
								m1v[2] = m2v[2];
								m2v[0] = tempAr[0];
								m2v[1] = tempAr[1];
								m2v[2] = tempAr[2];
							}
						}
					} // m2x closed
				} // m2y closed
			} // m2z closed
		}
	} //m1 closed
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL");

	//Finalize FFT

	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_FINALIZE");
	// FFT Finalize
	for (m1 = loop_min; m1 < loop_max; m1++) {
		if (curLevel > 1 or !_doNTGlobal or curLevel < _stopLevel) {
			if (_mpCellGlobalTop[curLevel][m1].occ == 0)
				continue;
		}
		radius = _mpCellGlobalTop[curLevel][m1].local.getRadius();
		FFTAccelerableExpansion &target =
										static_cast<bhfmm::SHLocalParticle &>(_mpCellGlobalTop[curLevel][m1].local).getExpansion();
		_FFTAcceleration->FFT_finalize_Target(target, radius);
#ifdef ENABLE_MPI
		if(curLevel > 1 && _doNTGlobal && curLevel >= _stopLevel){
			//iterate over all accessible cells m2 (they are always occ=0 -> extra treatment)
			m1v[0] = m1 % mpCells;
			m1v[1] = (m1 / mpCells) % mpCells;
			m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
			if(_numCellsOnGlobalLevel[1] == 2 && m1v[1] % 2 == 1){ //don't iterate over tower twice
				continue;
			}
			for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
				// to get periodic image
				m22z = (mpCells + m2z) % mpCells;

				m2v[2] = m2z;
				for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
					// to get periodic image
					m22y = (mpCells + m2y) % mpCells;

					m2v[1] = m2y;
					for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
						// to get periodic image
						m22x = (mpCells + m2x) % mpCells;
						if(_doNTGlobal and not(m2x == m1v[0] and m2z == m1v[2]) //tower
							and not(((_numCellsOnGlobalLevel[1] == 2 and m2y <= HiLim(1) - 2 and m2y >= LoLim(1) + 2) or m1v[1] == m2y)  and m2x >= LoLim(0) + 2 and not(m2x <= HiLim(0) - 2 and m2z < LoLim(2) + 2))){ //plate
							continue;
						}
						m2v[0] = m2x;

						m2 = (m22z * mpCells + m22y) * mpCells + m22x;

						if(_doNTGlobal and curLevel == 2 and ( m2y < LoLim(1) + 2 or m1 == m2 or (m2x > HiLim(0) - 2 and m2z < LoLim(2) + 2))){
							continue;
						}
						if(m1 == m2){
							continue;
						}
						radius = _mpCellGlobalTop[curLevel][m2].local.getRadius();
						FFTAccelerableExpansion& target3 =
								static_cast<bhfmm::SHLocalParticle&>(_mpCellGlobalTop[curLevel][m2].local).getExpansion();
						_FFTAcceleration->FFT_finalize_Target(target3, radius);

					}
				}
			}
		}
#endif
	}
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_FINALIZE");
//		std::cout << "global " <<n <<" all " << m <<" \n";
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL");
} // GatherWellSepLo_FFT_template closed

void UniformPseudoParticleContainer::M2LCompleteCell(int targetId, int level, int cellsPerDimension) {

	// corresponds to base_unit from UniformPseudoParticleContainer::M2LTowerPlateStep
	const double cellSize = 2.0 / sqrt(3);

	// get and initialize target
	int targetCoords[3];
	int cellsPerDimension2 = cellsPerDimension * cellsPerDimension;
	targetCoords[2] = targetId / cellsPerDimension2;
	targetCoords[1] = (targetId - (targetCoords[2] * cellsPerDimension2)) / cellsPerDimension;
	targetCoords[0] = targetId - (targetCoords[2] * cellsPerDimension2) - (targetCoords[1] * cellsPerDimension);

	FFTAccelerableExpansion& target = static_cast<bhfmm::SHLocalParticle&>(_mpCellGlobalTop[level][targetId]
			.local)
			.getExpansion();

	_FFTAcceleration->FFT_initialize_Target(target);

	// abort if cell target is empty
	if (_mpCellGlobalTop[level][targetId].occ == 0) {
		// Finalize
		double targetRadius = _mpCellGlobalTop[level][targetId].local.getRadius();
		_FFTAcceleration->FFT_finalize_Target(target, targetRadius);
		return;
	}

	int sourceCoords[3];
	// cell from where sources will be calculated
	int sourceRootCoords[3];
	// find start corner for stencil
	for (int i = 0; i < 3; ++i) {
		if (IsOdd(targetCoords[i]))
			sourceRootCoords[i] = targetCoords[i] - 3;
		else
			sourceRootCoords[i] = targetCoords[i] - 2;
	}

	// iterate through 6x6x6 grid around target cell with periodic boundaries
	for (int x = 0; x < 6; ++x) {
		sourceCoords[0] = sourceRootCoords[0] + x;
		if (sourceCoords[0] < 0)
			sourceCoords[0] += cellsPerDimension;
		else if (sourceCoords[0] >= cellsPerDimension)
			sourceCoords[0] -= cellsPerDimension;
		for (int y = 0; y < 6; ++y) {
			sourceCoords[1] = sourceRootCoords[1] + y;
			if (sourceCoords[1] < 0)
				sourceCoords[1] += cellsPerDimension;
			else if (sourceCoords[1] >= cellsPerDimension)
				sourceCoords[1] -= cellsPerDimension;
			for (int z = 0; z < 6; ++z) {
				sourceCoords[2] = sourceRootCoords[2] + z;
				if (sourceCoords[2] < 0)
					sourceCoords[2] += cellsPerDimension;
				else if (sourceCoords[2] >= cellsPerDimension)
					sourceCoords[2] -= cellsPerDimension;
				// skip cells that are too near
				if (std::abs(sourceRootCoords[0] + x - targetCoords[0]) < 2
					&& std::abs(sourceRootCoords[1] + y - targetCoords[1]) < 2
					&& std::abs(sourceRootCoords[2] + z - targetCoords[2]) < 2)
					continue;

				int sourceId = (sourceCoords[2] * cellsPerDimension + sourceCoords[1]) * cellsPerDimension
							   + sourceCoords[0];

				if(_mpCellGlobalTop[level][sourceId].occ == 0)
					continue;

				FFTAccelerableExpansion &source = static_cast<bhfmm::SHMultipoleParticle &>(_mpCellGlobalTop[level][sourceId]
						.multipole)
						.getExpansion();

				auto transferFunction = _FFT_TM->getTransferFunction(sourceRootCoords[0] + x - targetCoords[0],
																	 sourceRootCoords[1] + y - targetCoords[1],
																	 sourceRootCoords[2] + z - targetCoords[2],
																	 cellSize,
																	 cellSize,
																	 cellSize);

//                cout << "Level: " << level
//                     << " | Target: " << std::setw(5) << targetId
//                     << " (" << std::setw(2) << targetCoords[0]
//                     << ", " << std::setw(2) << targetCoords[1]
//                     << ", " << std::setw(2) << targetCoords[2]
//                     << ") | Source: " << std::setw(5) << sourceId
//                     << " (" << std::setw(2) << sourceCoords[0]
//                     << ", " << std::setw(2) << sourceCoords[1]
//                     << ", " << std::setw(2) << sourceCoords[2]
//                     << ") | tf: "
//                     << "(" << std::setw(2) << sourceRootCoords[0] + x - targetCoords[0]
//                     << ", " << std::setw(2) << sourceRootCoords[1] + y - targetCoords[1]
//                     << ", " << std::setw(2) << sourceRootCoords[2] + z - targetCoords[2]
//                     << ")" << std::endl;

				// calculate single M2L
				if(FFTSettings::USE_ORDER_REDUCTION){
					auto m2l_order = FFTOrderReduction::getM2LOrder(
							sourceRootCoords[2] + z - targetCoords[2],
							sourceRootCoords[1] + y - targetCoords[1],
							sourceRootCoords[0] + x - targetCoords[0],
							_maxOrd);
					static_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_OrderReduction_vec(
							source, target, transferFunction, m2l_order);

				} else
					_FFTAcceleration->FFT_M2L_vec(source, target, transferFunction);
			}
		}
	}

	// Finalize
	double targetRadius = _mpCellGlobalTop[level][targetId].local.getRadius();
	_FFTAcceleration->FFT_finalize_Target(target, targetRadius);
}

void UniformPseudoParticleContainer::M2LPair2Way(int cellA, int cellB, int level, int cellsPerDimension) {
	const double cellSize = 2.0 / sqrt(3);

	int  coordsA[3],
		 coordsB[3];
	int cellsPerDimension2 = cellsPerDimension * cellsPerDimension;

	coordsA[0] = cellA / cellsPerDimension2;
	coordsA[1] = (cellA - (coordsA[0] * cellsPerDimension2)) / cellsPerDimension;
	coordsA[2] = cellA - (coordsA[0] * cellsPerDimension2) - coordsA[1] * cellsPerDimension;

	coordsB[0] = cellB / cellsPerDimension2;
	coordsB[1] = (cellB - (coordsB[0] * cellsPerDimension2)) / cellsPerDimension;
	coordsB[2] = cellB - (coordsB[0] * cellsPerDimension2) - coordsB[1] * cellsPerDimension;

	FFTAccelerableExpansion &targetA = static_cast<bhfmm::SHLocalParticle &>(_mpCellGlobalTop[level][cellA]
			.local)
			.getExpansion();
	FFTAccelerableExpansion &sourceA = static_cast<bhfmm::SHMultipoleParticle &>(_mpCellGlobalTop[level][cellA]
			.multipole)
			.getExpansion();
	FFTAccelerableExpansion &targetB = static_cast<bhfmm::SHLocalParticle &>(_mpCellGlobalTop[level][cellB]
			.local)
			.getExpansion();
	FFTAccelerableExpansion &sourceB = static_cast<bhfmm::SHMultipoleParticle &>(_mpCellGlobalTop[level][cellB]
			.multipole)
			.getExpansion();

	auto transferFunction1 = _FFT_TM->getTransferFunction(coordsB[0] - coordsA[0],
														 coordsB[1] - coordsA[1],
														 coordsB[2] - coordsA[2],
														 cellSize,
														 cellSize,
														 cellSize);
	auto transferFunction2 = _FFT_TM->getTransferFunction(coordsA[0] - coordsB[0],
														 coordsA[1] - coordsB[1],
														 coordsA[2] - coordsB[2],
														 cellSize,
														 cellSize,
														 cellSize);

	// calculate single M2L
	if(FFTSettings::USE_ORDER_REDUCTION){
		auto m2l_order = FFTOrderReduction::getM2LOrder(
				coordsA[0] - coordsB[0],
				coordsA[1] - coordsB[1],
				coordsA[2] - coordsB[2],
				_maxOrd);
		// FIXME: two FFT_M2Ls are faster than one 2way
//		dynamic_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_2way_ORed_vec(
//				sourceA, sourceB, targetA, targetB, transferFunction1, m2l_order);
		dynamic_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_OrderReduction_vec(
				sourceA, targetB, transferFunction1, m2l_order);
		dynamic_cast<FFTAccelerationAPI_full*>(_FFTAcceleration)->FFT_M2L_OrderReduction_vec(
				sourceB, targetA, transferFunction2, m2l_order);

	} else {
//		dynamic_cast<FFTAccelerationAPI_2Way*>(_FFTAcceleration)->FFT_M2L_2way_vec(
//				sourceA, sourceB, targetA, targetB, transferFunction1);
		_FFTAcceleration->FFT_M2L_vec(sourceA, targetB, transferFunction1);
		_FFTAcceleration->FFT_M2L_vec(sourceB, targetA, transferFunction2);
	}
}
template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
void UniformPseudoParticleContainer::GatherWellSepLo_FFT_Local_template(
		double *cellWid, Vector3<int> localMpCells, int curLevel, int doHalos) {
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_LOKAL");
	//adjust for local level
	curLevel                  = curLevel - _globalLevel - 1;
	int             m1v[3];
	int             m2v[3];
	int             m1x, m1y, m1z;
	int             m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
	int             m22x, m22y, m22z; // for periodic image
	Vector3<double> periodicShift(0.0);
	//FFT param
	double          radius;
	double          base_unit = 2.0 / sqrt(3);
	int             M2L_order;
	FFTDataContainer *tf;
	//exclude halo values
	int zStart = 2;
	int xStart = 2;
	int yStart = 2;
	int xEnd   = localMpCells[0] - 2 - xStart;
	int yEnd   = localMpCells[1] - 2 - yStart;
	int zEnd   = localMpCells[2] - 2 - zStart;
	if (doHalos) { //in this case also iterate over halo values
		xStart = yStart = zStart = 0;
		xEnd   = localMpCells[0];
		yEnd   = localMpCells[1];
		zEnd   = localMpCells[2];
	}
	//Initialize FFT
	for (int mloop = 0; mloop < xEnd * yEnd * zEnd; mloop++) {
		m1x = mloop % xEnd + xStart;
		m1y = (mloop / xEnd) % yEnd + yStart;
		m1z = mloop / (xEnd * yEnd) + zStart;
		if (doHalos and not(m1z < 2 or m1z >= localMpCells[2] - 2) and not(m1y < 2 or m1y >= localMpCells[1] - 2) and
			not(m1x < 2 or m1x >= localMpCells[0] - 2)) {
			continue; //only initialize halo values in doHalos run
		}
		m1 = ((m1z) * localMpCells[1] + m1y) * localMpCells[0] + m1x;
		if (_mpCellLocal[curLevel][m1].occ == 0) { //initialize only cells that have occ=1
			continue;
		}
		radius = _mpCellLocal[curLevel][m1].local.getRadius();
		FFTAccelerableExpansion &source =
										static_cast<bhfmm::SHMultipoleParticle &>(_mpCellLocal[curLevel][m1].multipole).getExpansion();
		FFTAccelerableExpansion &target =
										static_cast<bhfmm::SHLocalParticle &>(_mpCellLocal[curLevel][m1].local).getExpansion();
		_FFTAcceleration->FFT_initialize_Source(source, radius);
		_FFTAcceleration->FFT_initialize_Target(target);
	}
	//calculate number of cells that need to be traversed
	int      numInnerCells[3];
	for (int d     = 0; d < 3; d++) {
		numInnerCells[d] = localMpCells[d] - 4;
	}
	int      numCellsToIterate[3];
	for (int d     = 0; d < 3; d++) {
		if (not(_doNTLocal and doHalos)) {
			numCellsToIterate[d] = numInnerCells[d];
		} else {
			numCellsToIterate[d] = localMpCells[d];
		}
	}
	//M2L in Fourier space
	for (int mloop = 0; mloop < numCellsToIterate[0] * numCellsToIterate[1] * numCellsToIterate[2]; mloop++) {
		int offset;
		if (_doNTLocal and doHalos) { //in NT method writing to halo is possible
			offset = 0;
		} else {
			offset = 2;
		}
		m1x = mloop % numCellsToIterate[0] + offset;
		m1y = (mloop / numCellsToIterate[0]) % numCellsToIterate[1] + offset;
		m1z = mloop / (numCellsToIterate[0] * numCellsToIterate[1]) + offset;
		m1 = ((m1z) * localMpCells[1] + m1y) * localMpCells[0] + m1x;
		if (filterM1Local(doHalos, m1, m1x, m1y, m1z, localMpCells, curLevel)) {
			continue;
		}
		m1v[0] = m1x;
		m1v[1] = m1y;
		m1v[2] = m1z;
		bool inHaloz, inHaloy, inHalox; //shows if current cell is in halo area in respective coordinate axis
		inHaloz  = inHaloy = inHalox = 0;
		for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
			if ((m2z < 2 or m2z >= localMpCells[2] - 2)) { //halo cell
				if (!doHalos) {
					continue;
				}
				if (doHalos) {
					inHaloz = 1;
				}
			}
			m2v[2] = m2z;
			for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
				if ((m2y < 2 or m2y >= localMpCells[1] - 2)) { //halo cell
					if (!doHalos) {
						continue;
					}
					if (doHalos) {
						if (_doNTLocal) {//m2 not in y halo allowed in NT (m2 in plate)
							continue;
						} else {
							inHaloy = 1;
						}
					}
				}
				m2v[1] = m2y;
				for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
					if ((m2x < 2 or m2x >= localMpCells[0] - 2)) { //halo cell
						if (!doHalos) {
							continue;
						}
						if (doHalos) {
							inHalox = 1;
						}
					}
					m2v[0] = m2x;
					m2      = (m2z * localMpCells[1] + m2y) * localMpCells[0] + m2x;
					if (filterM2Local(doHalos, m1, m1x, m1y, m1z, m2, m2x, m2y, m2z, localMpCells, curLevel, inHaloz,
									  inHaloy, inHalox)) {
						inHalox = 0;
						continue;
					}
					inHalox = 0;

					if (UseM2L_2way) {
						if (m1 > m2 && not(doHalos))
							continue;
					}
					int doBothDirections = (_doNTLocal and doHalos and not UseM2L_2way) ? 1
																						: 0; //do we need to calculate forces for both directions?

					for (int i = 0; i <= doBothDirections; i++) {

						tf = _FFT_TM->getTransferFunction(m2v[0] - m1v[0],
														  m2v[1] - m1v[1], m2v[2] - m1v[2], base_unit, base_unit,
														  base_unit);

						if (UseM2L_2way) {
							FFTAccelerableExpansion &source2 =
															static_cast<bhfmm::SHMultipoleParticle &>(_mpCellLocal[curLevel][m2].multipole).getExpansion();
							FFTAccelerableExpansion &target2 =
															static_cast<bhfmm::SHLocalParticle &>(_mpCellLocal[curLevel][m2].local).getExpansion();
							FFTAccelerableExpansion &source1 =
															static_cast<bhfmm::SHMultipoleParticle &>(_mpCellLocal[curLevel][m1].multipole).getExpansion();
							FFTAccelerableExpansion &target1 =
															static_cast<bhfmm::SHLocalParticle &>(_mpCellLocal[curLevel][m1].local).getExpansion();

							if (UseOrderReduction) {
								M2L_order = FFTOrderReduction::getM2LOrder(
										m2v[0] - m1v[0], m2v[1] - m1v[1],
										m2v[2] - m1v[2], _maxOrd);
								if (UseVectorization) {
									static_cast<FFTAccelerationAPI_full *>(_FFTAcceleration)->FFT_M2L_2way_ORed_vec(
											source2, source1, target2, target1, tf,
											M2L_order);
								} else {
									static_cast<FFTAccelerationAPI_full *>(_FFTAcceleration)->FFT_M2L_2way_ORed(
											source2, source1, target2, target1, tf,
											M2L_order);
								}
							} else {
								if (UseVectorization) {
									static_cast<FFTAccelerationAPI_2Way *>(_FFTAcceleration)->FFT_M2L_2way_vec(
											source2, source1, target2, target1, tf);
								} else {
									static_cast<FFTAccelerationAPI_2Way *>(_FFTAcceleration)->FFT_M2L_2way(
											source2, source1, target2, target1, tf);
								}
							}
						} else {
							FFTAccelerableExpansion &source =
															static_cast<bhfmm::SHMultipoleParticle &>(_mpCellLocal[curLevel][m2].multipole).getExpansion();
							FFTAccelerableExpansion &target =
															static_cast<bhfmm::SHLocalParticle &>(_mpCellLocal[curLevel][m1].local).getExpansion();
							if (UseOrderReduction) {
								M2L_order = FFTOrderReduction::getM2LOrder(
										m2v[0] - m1v[0], m2v[1] - m1v[1],
										m2v[2] - m1v[2], _maxOrd);
								if (UseVectorization) {
									static_cast<FFTAccelerationAPI_full *>(_FFTAcceleration)->FFT_M2L_OrderReduction_vec(
											source, target, tf, M2L_order);
								} else {
									static_cast<FFTAccelerationAPI_full *>(_FFTAcceleration)->FFT_M2L_OrderReduction(
											source, target, tf, M2L_order);
								}
							} else {
								if (UseVectorization) {
									_FFTAcceleration->FFT_M2L_vec(source, target,
																  tf);
								} else {
									_FFTAcceleration->FFT_M2L(source, target, tf);
								}
							}
						}

						if (!UseTFMemoization) {
							delete tf; //free useless memory
						}

						if (_doNTLocal and doHalos and not UseM2L_2way) {
							//exchange m1 and m2 for second direction
							int temp = m1;
							m1 = m2;
							m2 = temp;
							int tempAr[3];
							tempAr[0] = m1v[0];
							tempAr[1] = m1v[1];
							tempAr[2] = m1v[2];
							m1v[0]    = m2v[0];
							m1v[1]    = m2v[1];
							m1v[2]    = m2v[2];
							m2v[0]    = tempAr[0];
							m2v[1]    = tempAr[1];
							m2v[2]    = tempAr[2];
						}
					}
				} // m2x closed
				inHaloy = 0;
			} // m2y closed
			inHaloz  = 0;
		} // m2z closed

	} //mloop closed

	//Finalize FFT
	if (doHalos) {
		for (int mloop = 0; mloop < numCellsToIterate[0] * numCellsToIterate[1] * numCellsToIterate[2]; mloop++) {
			int offset;
			if (_doNTLocal) { //in NT method forces for halo cells have been computed! -> have to be finalized
				offset = 0;
			} else {
				offset = 2;
			}
			m1x = mloop % numCellsToIterate[0] + offset;
			m1y = (mloop / numCellsToIterate[0]) % numCellsToIterate[1] + offset;
			m1z = mloop / (numCellsToIterate[0] * numCellsToIterate[1]) + offset;

			m1 = ((m1z) * localMpCells[1] + m1y) * localMpCells[0] + m1x;
			if (_mpCellLocal[curLevel][m1].occ == 0) { // only finalize non empty cells
				continue;
			}

			radius = _mpCellLocal[curLevel][m1].local.getRadius();
			FFTAccelerableExpansion &target =
											static_cast<bhfmm::SHLocalParticle &>(_mpCellLocal[curLevel][m1].local).getExpansion();
			_FFTAcceleration->FFT_finalize_Target(target, radius);
		}
	}
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_LOKAL");
} // GatherWellSepLo_FFT_MPI_template closed

FFTAccelerationAPI *UniformPseudoParticleContainer::getFFTAcceleration() {
	return _FFTAcceleration;
}
#endif /* FMM_FFT */

void UniformPseudoParticleContainer::PropagateCellLo_Global(double */*cellWid*/, int mpCells, int curLevel){
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_GLOBAL");
	int m1v[3];
	int m2v[3];
	int iDir, m1, m1x, m1y, m1z, m2;
	int mpCellsN = 2 * mpCells;
	int loop_min = 0;
	int loop_max = mpCells * mpCells * mpCells;
	for (m1 = loop_min; m1 < loop_max; m1++) { //iterate over all global cells

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		m1x = m1v[0];
		m1y = m1v[1];
		m1z = m1v[2];

			//only do M2M for cells that are on the way up the tree from own cell on global tree;
			// other occ values are still 0 in global tree
			if (_mpCellGlobalTop[curLevel][m1].occ == 0){
				continue;
			}

		//iterate over 8 children of m1
		for (iDir = 0; iDir < 8; ++iDir) {
			m2v[0] = 2 * m1x;
			m2v[1] = 2 * m1y;
			m2v[2] = 2 * m1z;

			if (IsOdd(iDir)) m2v[0]     = m2v[0] + 1;
			if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
			if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

			m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];

			_mpCellGlobalTop[curLevel][m1].local.actOnLocalParticle(
					_mpCellGlobalTop[curLevel + 1][m2].local); //L2L operation
		} // iDir
	}
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_GLOBAL");
} // PropogateCellLo

void UniformPseudoParticleContainer::M2MCompleteCell(int targetId, int level, int cellsPerDim) {
	int sourceId,
		sourceCoords[3],
		targetCoords[3],
		cellsPerDim2 = cellsPerDim * cellsPerDim;
	targetCoords[2] = targetId / cellsPerDim2;
	targetCoords[1] = (targetId - (targetCoords[2] * cellsPerDim2)) / cellsPerDim;
	targetCoords[0] = targetId - (targetCoords[2] * cellsPerDim2) - (targetCoords[1] * cellsPerDim);

	auto *target = &_mpCellGlobalTop[level][targetId];
	for (int i = 0; i < 8; ++i) {
		sourceCoords[0] = targetCoords[0] * 2;
		sourceCoords[1] = targetCoords[1] * 2;
		sourceCoords[2] = targetCoords[2] * 2;

		if (IsOdd(i    )) ++sourceCoords[0];
		if (IsOdd(i / 2)) ++sourceCoords[1];
		if (IsOdd(i / 4)) ++sourceCoords[2];

		sourceId = (sourceCoords[2] * cellsPerDim * 2 + sourceCoords[1]) * cellsPerDim * 2
				   + sourceCoords[0];

		auto *source = &_mpCellGlobalTop[level + 1][sourceId];

		if(source->occ == 0)
			continue;

		//M2M operation
		target->occ += source->occ;
		target->multipole.addMultipoleParticle(source->multipole);
	}
}

void UniformPseudoParticleContainer::L2LCompleteCell(int sourceId, int level, int cellsPerDimension) {

	if (_mpCellGlobalTop[level][sourceId].occ == 0){
		return;
	}

	int targetId,
		sourceCoords[3],
		targetCoords[3],
		currentEdgeLength2 = cellsPerDimension * cellsPerDimension;
	sourceCoords[2] = sourceId / currentEdgeLength2;
	sourceCoords[1] = (sourceId - (sourceCoords[2] * currentEdgeLength2)) / cellsPerDimension;
	sourceCoords[0] = sourceId - (sourceCoords[2] * currentEdgeLength2) - (sourceCoords[1] * cellsPerDimension);


	for (int i = 0; i < 8; ++i) {
		targetCoords[0] = 2 * sourceCoords[0];
		targetCoords[1] = 2 * sourceCoords[1];
		targetCoords[2] = 2 * sourceCoords[2];

		if (IsOdd(i    )) targetCoords[0] = targetCoords[0] + 1;
		if (IsOdd(i / 2)) targetCoords[1] = targetCoords[1] + 1;
		if (IsOdd(i / 4)) targetCoords[2] = targetCoords[2] + 1;

		targetId = (targetCoords[2] * cellsPerDimension * 2 + targetCoords[1]) * cellsPerDimension * 2
				   + targetCoords[0];

		// L2L operation
		_mpCellGlobalTop[level][sourceId].local.actOnLocalParticle(
				_mpCellGlobalTop[level + 1][targetId].local);
	}
}

void UniformPseudoParticleContainer::L2PCompleteCell(int targetId) {
	processFarField(_leafContainer->getCells()[targetId]);
}

void UniformPseudoParticleContainer::P2MCompleteCell(int sourceId) {
	processMultipole(_leafContainer->getCells()[sourceId]);
}

void UniformPseudoParticleContainer::PropagateCellLo_Local(double* /*cellWid*/, Vector3<int> localMpCells, int curLevel, Vector3<int> offset){
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_LOKAL");
	//	int m1v[3];
	int m2v[3];

	int iDir, m1, m1x, m1y, m1z, m2;

	Vector3<int> localMpCellsN;
	for (int     i = 0; i < 3; i++) {
		localMpCellsN[i] = 2 * (localMpCells[i] - 4) + 4;
	}
	//correct length in case of globalLevel is reached
	Vector3<int> localMpCellsRow;
	std::vector<std::vector<MpCell> > *mpCellCurLevel;
	int curLevelp1 = curLevel - _globalLevel;

	if (curLevel == _globalLevel) {
		for (int i = 0; i < 3; i++) {
			localMpCellsRow[i] = (localMpCells[i] - 4) * _numProcessorsPerDim[i];
		}
		mpCellCurLevel = &_mpCellGlobalTop;

	} else {
		localMpCellsRow = localMpCells;
		mpCellCurLevel  = &_mpCellLocal;
		//adjust level to local tree
		curLevel        = curLevel - _globalLevel - 1;
	}
	int      numInnerCells[3];
	for (int d     = 0; d < 3; d++) {
		numInnerCells[d] = localMpCells[d] - 4;
	}
	for (int mloop = 0; mloop < numInnerCells[0] * numInnerCells[1] * numInnerCells[2]; mloop++) {
		m1x = mloop % numInnerCells[0];
		m1y = (mloop / numInnerCells[0]) % numInnerCells[1];
		m1z = mloop / (numInnerCells[0] * numInnerCells[1]);
		m1 = ((m1z + offset[2]) * localMpCellsRow[1] + m1y + offset[1]) * localMpCellsRow[0] + m1x + offset[0];

		if ((*mpCellCurLevel)[curLevel][m1].occ == 0) { //only iterate over non empty cells
			continue;
		}

		for (iDir = 0; iDir < 8; iDir++) { //iterate over 8 children of m1
			//adjust for halo
			m2v[0] = 2 * m1x + 2;
			m2v[1] = 2 * m1y + 2;
			m2v[2] = 2 * m1z + 2;

			if (IsOdd(iDir)) m2v[0]     = m2v[0] + 1;
			if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
			if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

			m2 = (m2v[2] * localMpCellsN[1] + m2v[1]) * localMpCellsN[0] + m2v[0];

			(*mpCellCurLevel)[curLevel][m1].local.actOnLocalParticle(
					_mpCellLocal[curLevelp1][m2].local); //L2L operation
		} // iDir
	}
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_LOKAL");
} // PropogateCellLo_MPI

void UniformPseudoParticleContainer::processMultipole(ParticleCellPointers& cell){
	int                               cellIndexV[3];
	std::vector<std::vector<MpCell> > *mpCellMaxLevel;
	int                               maxLevel;
	if (_maxLevel == _globalLevel) {
		mpCellMaxLevel = &_mpCellGlobalTop;
		maxLevel       = _maxLevel;
	} else {
		mpCellMaxLevel = &_mpCellLocal;
		maxLevel       = _maxLevel - _globalLevel - 1;
	}
	for (int i = 0; i < 3; i++) {
#if defined(ENABLE_MPI)
		if(_maxLevel == _globalLevel){
			cellIndexV[i] = rint(cell.getBoxMin(i) / _cellLength[i]);
		}
		else{
			cellIndexV[i] = rint((cell.getBoxMin(i)-_bBoxMin[i]) / _cellLength[i]) + 2;
			if(cellIndexV[i] < 2){
				std::cout << "Negative value " << cellIndexV[i] << " boxMin: " << cell.getBoxMin(i) << " bBoxMin: " << _bBoxMin[i] << " cellLength: " << _cellLength[i] << " \n";
			}
			if(cell.getBoxMin(i)-_bBoxMin[i]<0){
				std::cout << "Negative value\n";
			}
		}
#else
		cellIndexV[i] = rint(cell.getBoxMin(i) / _cellLength[i]);
#endif
	}
#if defined(ENABLE_MPI)
	int cellIndex;
	if (_maxLevel == _globalLevel){
		cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1]) * _globalNumCellsPerDim) + cellIndexV[0];
	}
	else{
		Vector3<int> numCellsOnLocalMaxLevel;
		for(int i = 0; i < 3; i++){
			numCellsOnLocalMaxLevel[i] = (_globalNumCellsPerDim / _numProcessorsPerDim[i] + 4);
		}
		cellIndex = ((numCellsOnLocalMaxLevel[1] * cellIndexV[2] + cellIndexV[1])	* numCellsOnLocalMaxLevel[0]) + cellIndexV[0];
	}
#else
	int cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1]) * _globalNumCellsPerDim) + cellIndexV[0];
#endif


//	mardyn_assert(cell.isInActiveWindow());

	int currentParticleCount = cell.getMoleculeCount();
	int Occupied             = 0;

	// loop over all particles in the cell
	for (int i = 0; i < currentParticleCount; i++) {
		++Occupied;
		Molecule &molecule1 = cell.moleculesAt(i);
		int      ni         = molecule1.numCharges();

		for (int j = 0; j < ni; j++) {
			const std::array<double, 3> dii      = molecule1.charge_d(j);
			const Charge                &chargei = static_cast<const Charge &> (molecule1.component()->charge(j));
			double                      dr[3];

			for (int k = 0; k < 3; k++) {
				dr[k] = molecule1.r(k) + dii[k];
			}    // for k closed

			bhfmm::Vector3<double> site_pos_vec3(dr);
			(*mpCellMaxLevel)[maxLevel][cellIndex].multipole.addSource(site_pos_vec3, chargei.q());

		}// for j closed
	} // current particle closed

	(*mpCellMaxLevel)[maxLevel][cellIndex].occ = Occupied;
}

void UniformPseudoParticleContainer::processFarField(ParticleCellPointers& cell) {
	int                               cellIndexV[3];
	std::vector<std::vector<MpCell> > *mpCellMaxLevel;

	int maxLevel;
	if (_maxLevel == _globalLevel) {
		mpCellMaxLevel = &_mpCellGlobalTop;
		maxLevel       = _maxLevel;
	} else {
		mpCellMaxLevel = &_mpCellLocal;
		maxLevel       = _maxLevel - _globalLevel - 1;
	}
	for (int i = 0; i < 3; i++) {
#if defined(ENABLE_MPI)
		if(_maxLevel == _globalLevel){
			cellIndexV[i] = rint(cell.getBoxMin(i) / _cellLength[i]);
		}
		else{
			cellIndexV[i] = rint((cell.getBoxMin(i)-_bBoxMin[i]) / _cellLength[i]) + 2;
		}
#else
		cellIndexV[i] = rint(cell.getBoxMin(i) / _cellLength[i]);
#endif
	}

#if defined(ENABLE_MPI)
	int cellIndex;
	if (_maxLevel == _globalLevel){
		cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1]) * _globalNumCellsPerDim) + cellIndexV[0];
	}
	else{
		Vector3<int> numCellsOnLocalMaxLevel;
		for(int i = 0; i < 3; i++){
			numCellsOnLocalMaxLevel[i] = (_globalNumCellsPerDim / _numProcessorsPerDim[i] + 4);
		}
		cellIndex = ((numCellsOnLocalMaxLevel[1] * cellIndexV[2] + cellIndexV[1])	* numCellsOnLocalMaxLevel[0]) + cellIndexV[0];
	}
#else
	int cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1]) * _globalNumCellsPerDim) + cellIndexV[0];
#endif

	SolidHarmonicsExpansion leLocal(_maxOrd);

	int             currentParticleCount = cell.getMoleculeCount();
	double          u                    = 0;
	double          uSum                 = 0.0;
	double          f[3]                 = {0.0, 0.0, 0.0};
	Vector3<double> f_vec3;
	double          virialSum            = 0.0;
	double          P_xxSum              = 0.0;
	double          P_yySum              = 0.0;
	double          P_zzSum              = 0.0;

	// loop over all particles in the cell
	for (int i = 0; i < currentParticleCount; i++) {
		Molecule &molecule1 = cell.moleculesAt(i);
		int      ni         = molecule1.numCharges();

		for (int j = 0; j < ni; j++) {
			const std::array<double, 3> dii      = molecule1.charge_d(j);
			const Charge                &chargei = static_cast<const Charge &> (molecule1.component()->charge(j));
			Vector3<double>             dr;

			for (int k = 0; k < 3; k++) {
				dr[k] = molecule1.r(k) + dii[k];
			}       // for k closed

			(*mpCellMaxLevel)[maxLevel][cellIndex].local.actOnTarget(dr, chargei.q(), u, f_vec3);
			f[0] = f_vec3[0];
			f[1] = f_vec3[1];
			f[2] = f_vec3[2];

			double   virial = 0.0;
			for (int l      = 0; l < 3; l++) {
				virial += -f[l] * dr[l];
			}
			P_xxSum += 0.5 * -f[0] * dr[0];
			P_yySum += 0.5 * -f[1] * dr[1];
			P_zzSum += 0.5 * -f[2] * dr[2];
			molecule1.Fchargeadd(j, f);
			uSum += 0.5 * u;
			virialSum += 0.5 * virial;
		}// for j closed
	} // current particle closed

		_domain->setLocalUpot(uSum + _domain->getLocalUpot());
		_domain->setLocalVirial(virialSum + _domain->getLocalVirial());
	//	_domain->addLocalP_xx(P_xxSum);
	//	_domain->addLocalP_yy(P_yySum);
	//	_domain->addLocalP_zz(P_zzSum);
}

void UniformPseudoParticleContainer::clear() {
	for (int n = 0; n < _maxLevel-_globalLevel; n++) {
		Vector3<int> localMpCells = _numCellsOnGlobalLevel * pow(2, n + 1) + Vector3<int>(4);

		for (int m1z = 0; m1z < localMpCells[2]; m1z++) {
			for (int m1y = 0; m1y < localMpCells[1]; m1y++) {
				for (int m1x = 0; m1x < localMpCells[0]; m1x++) {
					int cellIndexNew = (m1z * localMpCells[1] + m1y) * localMpCells[0] + m1x;

					_mpCellLocal[n][cellIndexNew].occ = 0;

					_mpCellLocal[n][cellIndexNew].multipole.clear();
					_mpCellLocal[n][cellIndexNew].local.clear();
				}
			}
		}
	}
	for (int n = _globalLevel; n >= 1; n--) {
		int mpCells = pow(2, n);

		for (int m1z = 0; m1z < mpCells; m1z++) {
			for (int m1y = 0; m1y < mpCells; m1y++) {
				for (int m1x = 0; m1x < mpCells; m1x++) {
					int cellIndexNew = (m1z * mpCells + m1y) * mpCells + m1x;

					_mpCellGlobalTop[n][cellIndexNew].occ = 0;

					_mpCellGlobalTop[n][cellIndexNew].multipole.clear();
					_mpCellGlobalTop[n][cellIndexNew].local.clear();
				}
			}
		}
	}

	// clear the MPI buffers
#ifdef ENABLE_MPI
	std::fill(_coeffVector.begin(), _coeffVector.end(), 0.0);
	std::fill(_occVector.begin(), _occVector.end(), 0);
#endif
}

void UniformPseudoParticleContainer::AllReduceMultipoleMoments() {
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE");
#ifdef ENABLE_MPI
	int coeffIndex = 0;
	for (int cellIndex = 0; cellIndex < _globalNumCells; cellIndex++) {
		const MpCell & currentCell = _mpCellGlobalTop[_maxLevel][cellIndex];

		// NOTE: coeffIndex modified in following call:
		currentCell.multipole.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

		mardyn_assert(cellIndex < _globalNumCells);
		_occVector[cellIndex] = currentCell.occ;
	}

	MPI_Allreduce(MPI_IN_PLACE, _coeffVector.data(), _coeffVectorLength*2, MPI_DOUBLE, MPI_SUM, _comm);
	MPI_Allreduce(MPI_IN_PLACE, _occVector.data(), _globalLevelNumCells, MPI_INT, MPI_SUM, _comm);


	coeffIndex = 0;
	for (int cellIndex = 0; cellIndex < _globalLevelNumCells; cellIndex++) {
		MpCell & currentCell = _mpCellGlobalTop[_maxLevel][cellIndex];

		currentCell.occ = _occVector[cellIndex];
		currentCell.multipole.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector.begin(), _coeffVector.end(), 0.0);

#endif
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE");
}

void UniformPseudoParticleContainer::AllReduceMultipoleMomentsLevelToTop(int numCellsLevel,int startingLevel) {
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE");
#ifdef ENABLE_MPI
	if(startingLevel > 0){
		int coeffIndex = 0;
		int numCellsLevelTemp = numCellsLevel;

		for(int level = startingLevel; level > 0 ; level--){
			for (int cellIndex = 0; cellIndex < numCellsLevelTemp; cellIndex++) {
				const MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];

				// NOTE: coeffIndex modified in following call:
				currentCell.multipole.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

			}
			numCellsLevelTemp /= 8;
		}
		int commLevel = (_avoidAllReduce && _stopLevel <= _globalLevel)? _stopLevel: _globalLevel;
		MPI_Iallreduce(MPI_IN_PLACE, _coeffVector.data(), _coeffVectorLength*2, MPI_DOUBLE, MPI_SUM, _allReduceComms[commLevel], &_allReduceRequest);
	}

#endif
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE");
}
void UniformPseudoParticleContainer::AllReduceMultipoleMomentsSetValues(int numCellsLevel,int startingLevel) {
#ifdef ENABLE_MPI
	if(startingLevel > 0){
	int coeffIndex = 0;
	int numCellsLevelTemp = numCellsLevel;
	numCellsLevelTemp = numCellsLevel;
	for(int level = startingLevel; level > 0 ; level--){
		for (int cellIndex = 0; cellIndex < numCellsLevelTemp; cellIndex++) {

			MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];
			currentCell.multipole.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

		}
		numCellsLevelTemp /= 8;
	}
	std::fill(_coeffVector.begin(), _coeffVector.end(), 0.0);
	}
#endif
}
void UniformPseudoParticleContainer::AllReduceLocalMoments(int mpCells, int _curLevel) {
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE_ME");

#ifdef ENABLE_MPI

	const int _row_Length=pow(mpCells, 3);
	int coeffIndex = 0;

	coeffIndex = 0;

	for (int cellIndex = 0; cellIndex < _row_Length; cellIndex++) {
		const MpCell & currentCell = _mpCellGlobalTop[_curLevel][cellIndex];

		if(currentCell.occ == 0) continue;
		currentCell.local.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

	}
	MPI_Allreduce(MPI_IN_PLACE, _coeffVector.data(), coeffIndex, MPI_DOUBLE, MPI_SUM, _comm);

	coeffIndex = 0;

	for (int cellIndex = 0; cellIndex < _row_Length; cellIndex++) {
		MpCell & currentCell = _mpCellGlobalTop[_curLevel][cellIndex];

		if(currentCell.occ == 0) continue;
		currentCell.local.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector.begin(), _coeffVector.end(), 0.0);

#endif
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE_ME");
}

void UniformPseudoParticleContainer::getHaloValues(Vector3<int> localMpCellsBottom,int bottomLevel, std::vector<double> buffer,
		int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh, bool doLocalExpansion){
#if defined(ENABLE_MPI)
	int coeffIndex = 0;
	Vector3<int> localMpCells = localMpCellsBottom;
	int cellIndex;
	int xLowLevel, yLowLevel, zLowLevel;
	int xHighLevel, yHighLevel, zHighLevel;
	for(int level=bottomLevel - _globalLevel - 1; level>= 0;level--){
		xLowLevel = (xLow < 0)? localMpCells[0] + xLow : xLow;
		yLowLevel = (yLow < 0)? localMpCells[1] + yLow : yLow;
		zLowLevel = (zLow < 0)? localMpCells[2] + zLow : zLow;
		xHighLevel = (xHigh <= 0)? localMpCells[0] + xHigh : xHigh;
		yHighLevel = (yHigh <= 0)? localMpCells[1] + yHigh : yHigh;
		zHighLevel = (zHigh <= 0)? localMpCells[2] + zHigh : zHigh;

		for (int z = zLowLevel; z < zHighLevel; z++) {
			for (int y = yLowLevel; y < yHighLevel; y++) {
				for (int x = xLowLevel; x < xHighLevel; x++) {
					cellIndex = (z * localMpCells[1] + y) * localMpCells[0] + x;
					const MpCell & currentCell = _mpCellLocal[level][cellIndex];
					if(doLocalExpansion){
						currentCell.local.writeValuesToMPIBuffer(buffer, coeffIndex);
					}
					else{
						currentCell.multipole.writeValuesToMPIBuffer(buffer, coeffIndex);
					}
				}
			}
		}
		for(int i = 0; i < 3; i++){
			localMpCells[i] = (localMpCells[i] - 4) / 2 + 4;
		}
	}
#endif
}

void UniformPseudoParticleContainer::setHaloValues(Vector3<int> localMpCellsBottom,int bottomLevel, std::vector<double> bufferRec,
		int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh, bool doLocalExpansion){
#if defined(ENABLE_MPI)

	int coeffIndex = 0;

	Vector3<int> localMpCells = localMpCellsBottom;
	int cellIndex;
	int xLowLevel, yLowLevel, zLowLevel;
	int xHighLevel, yHighLevel, zHighLevel;
	for(int level=bottomLevel - _globalLevel - 1; level>= 0;level--){
		xLowLevel = (xLow < 0)? localMpCells[0] + xLow : xLow;
		yLowLevel = (yLow < 0)? localMpCells[1] + yLow : yLow;
		zLowLevel = (zLow < 0)? localMpCells[2] + zLow : zLow;
		xHighLevel = (xHigh <= 0)? localMpCells[0] + xHigh : xHigh;
		yHighLevel = (yHigh <= 0)? localMpCells[1] + yHigh : yHigh;
		zHighLevel = (zHigh <= 0)? localMpCells[2] + zHigh : zHigh;
		for (int z = zLowLevel; z < zHighLevel; z++) {
			for (int y = yLowLevel; y < yHighLevel; y++) {
				for (int x = xLowLevel; x < xHighLevel; x++) {
					cellIndex = (z * localMpCells[1] + y) * localMpCells[0] + x;
					MpCell & currentCell = _mpCellLocal[level][cellIndex];
					if(doLocalExpansion){
						currentCell.local.addValuesFromMPIBuffer(bufferRec, coeffIndex);
					}
					else{
						currentCell.multipole.readValuesFromMPIBuffer(bufferRec, coeffIndex);
					}
					int empty = 1;
					for (int l = 0; l <= _maxOrd; ++l) {
						for (int m=0; m <= l; ++m){
							if(currentCell.multipole.getExpansion().getC(l,m) != 0.0){
								empty = 0;
								break;
							}
							if(currentCell.multipole.getExpansion().getS(l,m) != 0.0){
								empty = 0;
								break;
							}
						}
					}
					if(empty){
						currentCell.occ = 0;
					}
					else{
						currentCell.occ = 1;
					}

				}
			}
		}
		for(int i = 0; i < 3; i++){
			localMpCells[i] = (localMpCells[i] - 4) / 2 + 4;
		}
	}

#endif
}

void UniformPseudoParticleContainer::communicateHalos(){
	global_simulation->timers()->start("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMMUNICATION_HALOS");
	if(!_overlapComm){
		communicateHalosNoOverlap();
	}
	else{
		communicateHalosOverlapStart();
	}
	global_simulation->timers()->stop("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMMUNICATION_HALOS");
}
#ifdef ENABLE_MPI
void UniformPseudoParticleContainer::communicateOwnGlobalValue(int stopLevel, bool receive ){
	//set Values
	int bufferPosition = 0;
	//coords of MPI rank in cartesian grid
	int coords[3];
	int myRank;
	MPI_Comm_rank(_comm,&myRank);
	MPI_Cart_coords(_comm, myRank, 3, coords);
	int numRanks;
	MPI_Comm_size(_comm, &numRanks);
	if(_fuseGlobalCommunication){
		stopLevel = (stopLevel < 2)? 2 : stopLevel;
	}
	for(int level = _globalLevel; level >= stopLevel; level --){
		int mpCells = pow(2,level);
		int stride = pow(2,_globalLevel - level);
		//cell for which MPI rank is responsible on current level
		int coordsLevel[3];
		//min coords of rank in parent area
//		int coordsFloored[3];
		int coordsFlooredLevel[3];
		int cellIndex;
		for(int d = 0; d < 3; d++){
			coordsLevel[d] = ((coords[d] * _numCellsOnGlobalLevel[d]) / (stride));
			coordsFlooredLevel[d] = ((coordsLevel[d]/ 2) * 2);
		}
		int xIndex = coordsLevel[0];
		int yIndex = coordsLevel[1];
		int zIndex = coordsLevel[2];
		if(not receive){

			if(level == _globalLevel or _fuseGlobalCommunication){
				int endx, endy, endz;
				if(_fuseGlobalCommunication){
					endx = 2;
					endy = 2;
					endz = 2;
				}
				else{
					endx = _numCellsOnGlobalLevel[0];
					endy = _numCellsOnGlobalLevel[1];
					endz = _numCellsOnGlobalLevel[2];

				}
				int indexPosition = 0;
				for(int x = 0; x < endx; x ++){
					for(int y = 0; y < endy; y++){
						for(int z = 0; z < endz; z++){
							if(!_fuseGlobalCommunication){
								indexPosition = 0;
								cellIndex = ((zIndex + z) * mpCells + yIndex + y) * mpCells + xIndex + x;
							}
							else{
								cellIndex = ((coordsFlooredLevel[2] + z) * mpCells + coordsFlooredLevel[1] + y) * mpCells + coordsFlooredLevel[0] + x;
							}

							MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];
							if(!_fuseGlobalCommunication){
								currentCell.multipole.writeValuesToMPIBuffer(_multipoleBufferOverlapGlobal->getAreaBuffers()[8*(_globalLevel - level)+ z*4 + y*2 + x], indexPosition);
							}
							else{
								currentCell.multipole.writeValuesToMPIBuffer(_multipoleBufferOverlapGlobal->getAreaBuffers()[8*(_globalLevel - level)], indexPosition);
							}
						}
					}
				}
			}
			else{
				int indexPosition = 0;

				cellIndex = ((zIndex) * mpCells + yIndex) * mpCells + xIndex;
				MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];
				currentCell.multipole.writeValuesToMPIBuffer(_multipoleBufferOverlapGlobal->getAreaBuffers()[8*(_globalLevel - level)], indexPosition);

			}
		}

		else{ //can only happen with NT on global tree
			int start, end;
			if(_fuseGlobalCommunication){
				start = -1;
				end = 1;
			}
			else{
				start = -2;
				end = 3;
			}
			for(int x = start; x <= end; x++ ){
				for(int y = start; y <= end; y++){
					for(int z = start; z <= end; z++){
						//allow cells only if in plate or tower
						bool condition;
						if(!_fuseGlobalCommunication){
							if(level == _globalLevel){
								condition = ( floor((x + coordsFlooredLevel[0])/(1.0 * _numCellsOnGlobalLevel[0])) == coords[0] and floor((z + coordsFlooredLevel[2])/(1.0 * _numCellsOnGlobalLevel[2])) == coords[2]) //tower
										or (((y < 2 and y >= 0 and _numCellsOnGlobalLevel[1] == 2) or floor((y  + coordsFlooredLevel[1])/(1.0 * _numCellsOnGlobalLevel[1])) == coords[1]) and x < 2 and not (x >= 0 and z >= 2)); //plate (reversed for send)
							}
							else{
								condition = ( x + coordsFlooredLevel[0] == coordsLevel[0] and z + coordsFlooredLevel[2] == coordsLevel[2]) //tower
										or ( y  + coordsFlooredLevel[1] == coordsLevel[1] and x < 2 and not (x >= 0 and z >= 2)); //plate (reversed for send)

							}
						}
						else{
							condition = (x == 0 and z == 0) // tower
										or (y == 0 and x < 1 and not ( x==0 and z == 1)); //plate (reversed for send)
						}
//						condition = condition and (x<0 or x >= 2 or y < 0 or y >= 2 or z < 0 or z >= 2); //do not send within parent cell

						if(condition){
							bool filterDoubles;//filters out double accesses to procs
							if(!_fuseGlobalCommunication){
								filterDoubles = level > 2 or (level == 2 and y < 2 and not (x < 0 and z >= 2)) or (level == 1 and x < 2 and x >= 0 and y < 2 and y >= 0 and z < 2 and z >= 0); //filters out double accesses to procs
								int xTemp = (coordsFlooredLevel[0] + x + mpCells) % mpCells;
								int yTemp = (coordsFlooredLevel[1] + y + mpCells) % mpCells;
								int zTemp = (coordsFlooredLevel[2] + z + mpCells) % mpCells;
								filterDoubles = filterDoubles and not(xIndex == xTemp and yIndex == yTemp and zIndex == zTemp);
								//check if source of value is own rank -> needs to be filtered out
								if(level == _globalLevel){
									for(int i = 0; i < _numCellsOnGlobalLevel[0]; i++){
										for(int j = 0; j < _numCellsOnGlobalLevel[1]; j++){
											for(int k = 0; k < _numCellsOnGlobalLevel[2]; k++){
												filterDoubles = filterDoubles and not(xIndex == xTemp - i and yIndex == yTemp -j and zIndex == zTemp - k);
											}
										}

									}
								}

							}
							else{
								filterDoubles = level > 2 or (level == 2 and y < 1 and not (x == -1 and z == 1)); //filters out double accesses to procs
								filterDoubles = filterDoubles and not (x == 0 and y == 0 and z == 0);
							}


							if(filterDoubles){
								int indexPosition = 0;
								if(!_fuseGlobalCommunication or level != _globalLevel){
									//calculate offsets if more than one cell on global level
									const int xOffset = (level == _globalLevel) ? abs(x) % _numCellsOnGlobalLevel[0] : 0;
									const int yOffset = (level == _globalLevel) ? abs(y) % _numCellsOnGlobalLevel[1] : 0;
									const int zOffset = (level == _globalLevel) ? abs(z) % _numCellsOnGlobalLevel[2] : 0;
									cellIndex = ((zIndex + zOffset) * mpCells + yIndex + yOffset) * mpCells + xIndex + xOffset;
									MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];
									currentCell.local.addValuesFromMPIBuffer(_multipoleRecBufferOverlapGlobal->getAreaBuffers()[bufferPosition], indexPosition);
								}
								else{
									for(int xLocal = 0; xLocal < _numCellsOnGlobalLevel[0]; xLocal++){
										for(int yLocal = 0; yLocal < _numCellsOnGlobalLevel[1]; yLocal++){
											for(int zLocal = 0; zLocal < _numCellsOnGlobalLevel[2]; zLocal++){
												const int xIndexTemp = xIndex + xLocal;
												const int yIndexTemp = yIndex + yLocal;
												const int zIndexTemp = zIndex + zLocal;
												cellIndex = (zIndexTemp * mpCells + yIndexTemp) * mpCells + xIndexTemp;
												MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];
												currentCell.local.addValuesFromMPIBuffer(_multipoleRecBufferOverlapGlobal->getAreaBuffers()[bufferPosition], indexPosition);
											}
										}
									}
								}
							}
							bufferPosition++;
						}

					}
				}
			}
		}
	}
	//start communicating Values
	if(_globalLevel >= 1 and not receive and not (_globalLevel == 1 and _fuseGlobalCommunication)){
		_multipoleBufferOverlapGlobal->communicateGlobalLevels(_globalLevel,stopLevel);
	}
}

void UniformPseudoParticleContainer::communicateHaloGlobalValues(int stopLevel, bool send){
	//set Values
	int myRank;
	//coords of MPI rank in cartesian grid
	int coords[3];
	int bufferPosition = 0;
	MPI_Comm_rank(_comm,&myRank);
	MPI_Cart_coords(_comm, myRank, 3, coords);
	if(_fuseGlobalCommunication){
		stopLevel = (stopLevel < 2)? 2 : stopLevel;
	}
	for(int level = _globalLevel; level >= stopLevel; level --){
		int mpCells = pow(2,level);
		int stride = pow(2,_globalLevel - level);

		//min coords of MPI rank in parent area
		int coordsFloored[3];
		//min coords of cell in parent area
		int coordsFlooredLevel[3];
		//cell for which MPI rank is responsible on current level
		int coordsLevel[3];

		int cellIndex;
		for(int d = 0; d < 3; d++){
			coordsFloored[d] = ((coords[d] * _numCellsOnGlobalLevel[d]) / (2 * stride)) * 2 * stride;
			coordsLevel[d] = ((coords[d] * _numCellsOnGlobalLevel[d]) / (stride));
			coordsFlooredLevel[d] = ((coordsLevel[d]/ 2) * 2);


		}
		int start, end;
		if(_fuseGlobalCommunication){
			start = -1;
			end = 1;
		}
		else{
			start = -2;
			end = 3;
		}
		for(int x = start; x <= end; x++ ){
			for(int y = start; y <= end; y++){
				for(int z = start; z <= end; z++){
					int indexPosition = 0;
					bool condition;
					if(_doNTGlobal){ //allow cells only if in plate or tower
						if(!_fuseGlobalCommunication){
							if(level == _globalLevel){
								condition = ( floor((x + coordsFlooredLevel[0])/(1.0 * _numCellsOnGlobalLevel[0])) == coords[0] and floor((z + coordsFlooredLevel[2])/(1.0 * _numCellsOnGlobalLevel[2])) == coords[2]) //tower
										or (((y < 2 and y >= 0 and _numCellsOnGlobalLevel[1] == 2) or floor((y + coordsFlooredLevel[1])/(1.0 * _numCellsOnGlobalLevel[1])) == coords[1]) and x >= 0 and not (x < 2 and z < 0)); //plate

							}
							else{
								condition = ( x + coordsFlooredLevel[0] == coordsLevel[0] and z + coordsFlooredLevel[2] == coordsLevel[2]) //tower
																or (y + coordsFlooredLevel[1] == coordsLevel[1] and x >= 0 and not (x < 2 and z < 0)); //plate

							}
						}
						else{
							condition = (x == 0 and z == 0) // tower
										or (y == 0 and x >= 0 and not ( x < 1 and z == -1)); //plate
						}
					}
					else{
						if(!_fuseGlobalCommunication){
							condition = _importWholeGlobalRegion or abs(x +(coordsFloored[0] - coords[0])/stride) >= 2 or abs(y +(coordsFloored[1] - coords[1])/stride) >= 2 or abs(z + (coordsFloored[2] - coords[2])/stride) >= 2;
						}
						else{
							condition = x != 0 or y != 0 or z != 0;
						}
					}
					if(condition){
						if(!_fuseGlobalCommunication){
							int xIndex = (coordsFlooredLevel[0] + x + mpCells) % mpCells;
							int yIndex = (coordsFlooredLevel[1] + y + mpCells) % mpCells;
							int zIndex = (coordsFlooredLevel[2] + z + mpCells) % mpCells;
							cellIndex = (zIndex * mpCells + yIndex) * mpCells + xIndex;
							MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];
							if(send){
								currentCell.local.writeValuesToMPIBuffer(_multipoleBufferOverlapGlobal->getAreaBuffers()[bufferPosition], indexPosition);

							}
							else{
								currentCell.multipole.readValuesFromMPIBuffer(_multipoleRecBufferOverlapGlobal->getAreaBuffers()[bufferPosition], indexPosition);
							}
						}
						else{
							if(send){
								//offsets are 0 if more than one cell is on global level in this dimension is owned by processor
								int xOffset = coordsLevel[0] % 2;
								int yOffset = coordsLevel[1] % 2;
								int zOffset = coordsLevel[2] % 2;
								int endX = (level == _globalLevel)?  _numCellsOnGlobalLevel[0] : 1;
								int endY = (level == _globalLevel)?  _numCellsOnGlobalLevel[1] : 1;
								int endZ = (level == _globalLevel)?  _numCellsOnGlobalLevel[2] : 1;

								for(int xLocal = 0; xLocal < endX; xLocal++){
									for(int yLocal = 0; yLocal < endY; yLocal++){
										for(int zLocal = 0; zLocal < endZ; zLocal++){
											int xIndex = (coordsFlooredLevel[0] + 2 * x + xLocal + xOffset + mpCells) % mpCells;
											int yIndex = (coordsFlooredLevel[1] + 2 * y + yLocal + yOffset + mpCells) % mpCells;
											int zIndex = (coordsFlooredLevel[2] + 2 * z + zLocal + zOffset + mpCells) % mpCells;
											cellIndex = (zIndex * mpCells + yIndex) * mpCells + xIndex;
											MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];
											currentCell.local.writeValuesToMPIBuffer(_multipoleBufferOverlapGlobal->getAreaBuffers()[bufferPosition], indexPosition);
										}
									}
								}
							}
							else{
								for(int xLocal = 0; xLocal < 2; xLocal++){
									for(int yLocal = 0; yLocal < 2; yLocal++){
										for(int zLocal = 0; zLocal < 2; zLocal++){
											int xIndex = (coordsFlooredLevel[0] + 2 * x + xLocal+ mpCells) % mpCells;
											int yIndex = (coordsFlooredLevel[1] + 2 * y + yLocal + mpCells) % mpCells;
											int zIndex = (coordsFlooredLevel[2] + 2 * z + zLocal + mpCells) % mpCells;
											cellIndex = (zIndex * mpCells + yIndex) * mpCells + xIndex;
											MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];
											currentCell.multipole.readValuesFromMPIBuffer(_multipoleRecBufferOverlapGlobal->getAreaBuffers()[bufferPosition], indexPosition);
										}
									}
								}
							}
						}

						bufferPosition++;

					}
				}
			}
		}

	}
	if(_globalLevel >= 1 and send and _doNTGlobal and not (_globalLevel == 1 and _fuseGlobalCommunication)){ //backCommunication
		_multipoleBufferOverlapGlobal->communicateGlobalLevels(_globalLevel,stopLevel,true);
	}
}
#endif
void UniformPseudoParticleContainer::communicateHalosNoOverlap(){ //Outdated!!!

#if defined(ENABLE_MPI)
	_multipoleBuffer->clear();

	Vector3<int> localMpCellsBottom;
	for(int i = 0; i < 3; i++){
		localMpCellsBottom[i] = pow(2,_maxLevel) / _numProcessorsPerDim[i]  + 4;
	}
	//communicate along x axis
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getLeftBuffer(),
			2, 4, 2, -2, 2, -2, false);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getRightBuffer(),
			-4, -2, 2, -2, 2, -2, false);

	communicateHalosAlongAxis(_multipoleBuffer->getLeftBuffer(),_multipoleBuffer->getRightBuffer(),_multipoleRecBuffer->getLeftBuffer(),_multipoleRecBuffer->getRightBuffer(),
			_neighbours[0],_neighbours[1],_multipoleBuffer->getXSize()
			);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getLeftBuffer(),
			0, 2, 2, -2, 2, -2, false);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getRightBuffer(),
			-2, 0, 2, -2, 2, -2, false);

	//communicate along y axis
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getBottomBuffer(),
			0, 0, 2, 4, 2, -2, false);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getTopBuffer(),
			0, 0, -4, -2, 2, -2, false);

	communicateHalosAlongAxis(_multipoleBuffer->getBottomBuffer(),_multipoleBuffer->getTopBuffer(),_multipoleRecBuffer->getBottomBuffer(),_multipoleRecBuffer->getTopBuffer(),
			_neighbours[2],_neighbours[3],_multipoleBuffer->getYSize()
			);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getBottomBuffer(),
			0, 0, 0, 2, 2, -2, false);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getTopBuffer(),
			0, 0, -2, 0, 2, -2, false);

	//communicate along z axis
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getBackBuffer(),
			0, 0, 0, 0, 2, 4, false);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getFrontBuffer(),
			0, 0, 0, 0, -4, -2, false);

	communicateHalosAlongAxis(_multipoleBuffer->getBackBuffer(),_multipoleBuffer->getFrontBuffer(),_multipoleRecBuffer->getBackBuffer(),_multipoleRecBuffer->getFrontBuffer(),
			_neighbours[4],_neighbours[5],_multipoleBuffer->getZSize()
			);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getBackBuffer(),
			0, 0, 0, 0, 0, 2, false);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getFrontBuffer(),
			0, 0, 0, 0, -2, 0, false);


	MPI_Barrier(_comm);
#endif
}

void UniformPseudoParticleContainer::communicateHalosOverlapStart(){
#if defined(ENABLE_MPI)
	//clear buffers
	_multipoleBufferOverlap->clear();

	Vector3<int> localMpCellsBottom;
	for(int i = 0; i < 3; i++){
		localMpCellsBottom[i] = pow(2,_maxLevel) / _numProcessorsPerDim[i]  + 4;
	}
	//fill buffers of the halo areas of the simulation cube

	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[0],
										2, 4, 2, -2, 2, -2, false);
	if(not(_doNTLocal)){

		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[1],
									-4, -2, 2, -2, 2, -2, false);
	}
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[2],
									2, -2, 2, 4, 2, -2, false);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[3],
									2, -2, -4, -2, 2, -2, false);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[4],
										2, -2, 2, -2, 2, 4, false);
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[5],
									2, -2, 2, -2, -4, -2, false);
	}
	//fill edges buffers of the halo areas
	//adjacent edges to lower x area
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[0],
										2, 4, 2, 4, 2, -2, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[2],
										2, 4, -4, -2, 2, -2, false);
	}
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[4],
									2, 4, 2, -2, 2, 4, false);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[6],
									2, 4, 2, -2, -4, -2, false);


	//adjacent edges to higher x area
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[1],
										-4, -2, -4, -2, 2, -2, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[3],
										-4, -2, 2, 4, 2, -2, false);

		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[5],
										-4, -2, 2, -2, -4, -2, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[7],
									-4, -2, 2, -2, 2, 4, false);
	}
	//remaining edges
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[8],
										2, -2, 2, 4, 2, 4, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[10],
										2, -2, 2, 4, -4, -2, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[11],
										2, -2, -4, -2, 2, 4, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[9],
										2, -2, -4, -2, -4, -2, false);
	}
	//corners
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[0],
											2, 4, 2, 4, 2, 4, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[1],
											-4, -2, -4, -2, -4, -2, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[2],
											2, 4, 2, 4, -4, -2, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[3],
											-4, -2, -4, -2, 2, 4, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[4],
											2, 4, -4, -2, 2, 4, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[5],
											-4, -2, 2, 4, -4, -2, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[6],
											2, 4, -4, -2, -4, -2, false);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[7],
											-4, -2, 2, 4, 2, 4, false);
	}

	//start sending
//	_multipoleBufferOverlap->startCommunication();
	_multipoleBufferOverlap->communicate(false);
#endif
}

void UniformPseudoParticleContainer::communicateHalosOverlapSetHalos(){
#if defined(ENABLE_MPI)
	//read buffers

	Vector3<int> localMpCellsBottom;
	for(int i = 0; i < 3; i++){
		localMpCellsBottom[i] = pow(2,_maxLevel) / _numProcessorsPerDim[i]  + 4;
	}
	//read buffers of the halo areas of the simulation cube

	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[0],
									0, 2, 2, -2, 2, -2, false);
	}
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[1],
									-2, 0, 2, -2, 2, -2, false);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[2],
									2, -2, 0, 2, 2, -2, false);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[3],
									2, -2, -2, 0, 2, -2, false);
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[4],
									2, -2, 2, -2, 0, 2, false);
	}
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[5],
									2, -2, 2, -2, -2, 0, false);

	//read edges buffers of the halo areas
	//adjacent edges to lower x area
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[0],
										0, 2,0,2,2,-2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[2],
										0, 2,-2,0,2,-2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[4],
										0, 2, 2,-2, 0, 2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[6],
										0, 2,2,-2,-2, 0, false);
	}
	//adjacent edges to higher x area
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[1],
										-2, 0, -2, 0, 2,-2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[3],
										-2, 0, 0, 2, 2,-2, false);
	}
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[5],
									-2, 0, 2, -2, -2, 0, false);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[7],
									-2, 0, 2, -2, 0, 2, false);

	//remaining edges
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[8],
										2, -2, 0, 2, 0, 2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[10],
										2, -2, 0, 2, -2, 0, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[11],
										2, -2, -2, 0, 0, 2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[9],
										2, -2, -2, 0, -2, 0, false);
	}

	//corners
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[0],
											0, 2, 0, 2, 0, 2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[1],
											-2, 0, -2, 0, -2, 0, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[2],
											0, 2, 0, 2, -2, 0, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[3],
											-2, 0, -2, 0, 0, 2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[4],
											0, 2, -2, 0, 0, 2, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[5],
											-2, 0, 0, 2, -2, 0, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[6],
											0, 2, -2, 0, -2, 0, false);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[7],
											-2, 0, 0, 2, 0, 2, false);
	}
#endif
}

void UniformPseudoParticleContainer::communicateHalosOverlapPostProcessingSetHalos(){
#if defined(ENABLE_MPI)
	//start receiving



	Vector3<int> localMpCellsBottom;
	for(int i = 0; i < 3; i++){
		localMpCellsBottom[i] = pow(2,_maxLevel) / _numProcessorsPerDim[i]  + 4;
	}
	//fill buffers of the halo areas of the simulation cube

	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[0],
											2, 4, 2, -2, 2, -2, true);
	if(not(_doNTLocal)){

		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[1],
									-4, -2, 2, -2, 2, -2, true);
	}
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[2],
									2, -2, 2, 4, 2, -2, true);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[3],
									2, -2, -4, -2, 2, -2, true);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[4],
										2, -2, 2, -2, 2, 4, true);
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[5],
									2, -2, 2, -2, -4, -2, true);
	}
	//fill edges buffers of the halo areas
	//adjacent edges to lower x area
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[0],
										2, 4, 2, 4, 2, -2, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[2],
										2, 4, -4, -2, 2, -2, true);
	}
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[4],
									2, 4, 2, -2, 2, 4, true);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[6],
									2, 4, 2, -2, -4, -2, true);


	//adjacent edges to higher x area
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[1],
										-4, -2, -4, -2, 2, -2, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[3],
										-4, -2, 2, 4, 2, -2, true);

		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[5],
										-4, -2, 2, -2, -4, -2, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[7],
									-4, -2, 2, -2, 2, 4, true);
	}
	//remaining edges
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[8],
										2, -2, 2, 4, 2, 4, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[10],
										2, -2, 2, 4, -4, -2, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[11],
										2, -2, -4, -2, 2, 4, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[9],
										2, -2, -4, -2, -4, -2, true);
	}
	//corners
	if(not(_doNTLocal)){
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[0],
											2, 4, 2, 4, 2, 4, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[1],
											-4, -2, -4, -2, -4, -2, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[2],
											2, 4, 2, 4, -4, -2, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[3],
											-4, -2, -4, -2, 2, 4, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[4],
											2, 4, -4, -2, 2, 4, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[5],
											-4, -2, 2, 4, -4, -2, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[6],
											2, 4, -4, -2, -4, -2, true);
		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[7],
											-4, -2, 2, 4, 2, 4, true);
	}
#endif
}

void UniformPseudoParticleContainer::communicateHalosOverlapPostProcessingStart(){
#if defined(ENABLE_MPI)
	//clear buffers
	_multipoleBufferOverlap->clear();
	//read buffers

	Vector3<int> localMpCellsBottom;
	for(int i = 0; i < 3; i++){
		localMpCellsBottom[i] = pow(2,_maxLevel) / _numProcessorsPerDim[i]  + 4;
	}
	//read buffers of the halo areas of the simulation cube

	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[0],
									0, 2, 2, -2, 2, -2, true);
	}
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[1],
									-2, 0, 2, -2, 2, -2, true);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[2],
									2, -2, 0, 2, 2, -2, true);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[3],
									2, -2, -2, 0, 2, -2, true);
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[4],
									2, -2, 2, -2, 0, 2, true);
	}
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[5],
									2, -2, 2, -2, -2, 0, true);

	//read edges buffers of the halo areas
	//adjacent edges to lower x area
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[0],
										0, 2,0,2,2,-2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[2],
										0, 2,-2,0,2,-2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[4],
										0, 2, 2,-2, 0, 2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[6],
										0, 2,2,-2,-2, 0, true);
	}
	//adjacent edges to higher x area
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[1],
										-2, 0, -2, 0, 2,-2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[3],
										-2, 0, 0, 2, 2,-2, true);
	}
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[5],
									-2, 0, 2, -2, -2, 0, true);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[7],
									-2, 0, 2, -2, 0, 2, true);

	//remaining edges
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[8],
										2, -2, 0, 2, 0, 2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[10],
										2, -2, 0, 2, -2, 0, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[11],
										2, -2, -2, 0, 0, 2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[9],
										2, -2, -2, 0, -2, 0, true);
	}

	//corners
	if(not(_doNTLocal)){
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[0],
											0, 2, 0, 2, 0, 2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[1],
											-2, 0, -2, 0, -2, 0, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[2],
											0, 2, 0, 2, -2, 0, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[3],
											-2, 0, -2, 0, 0, 2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[4],
											0, 2, -2, 0, 0, 2, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[5],
											-2, 0, 0, 2, -2, 0, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[6],
											0, 2, -2, 0, -2, 0, true);
		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[7],
											-2, 0, 0, 2, 0, 2, true);
	}
	//start sending
//	_multipoleBufferOverlap->startCommunication();
	_multipoleBufferOverlap->communicate(true);
#endif
}

void UniformPseudoParticleContainer::communicateHalosAlongAxis(std::vector<double>lowerNeighbourBuffer,
															   std::vector<double>higherNeighbourBuffer,
															   std::vector<double>lowerNeighbourBufferRec,
															   std::vector<double>higherNeighbourBufferRec,
															   int lowerNeighbour,
															   int higherNeighbour,
															   int haloSize) {
#if defined(ENABLE_MPI)
	MPI_Request low, high;

	MPI_Status lowRecv,highRecv;

	MPI_Isend(lowerNeighbourBuffer.data(), haloSize, MPI_DOUBLE, lowerNeighbour, 1,
	_comm, &low);
	MPI_Isend(higherNeighbourBuffer.data(), haloSize, MPI_DOUBLE, higherNeighbour, 3,
	_comm, &high);

	MPI_Recv(lowerNeighbourBufferRec.data(), haloSize,MPI_DOUBLE, lowerNeighbour,3,_comm, &lowRecv);
	MPI_Recv(higherNeighbourBufferRec.data(), haloSize,MPI_DOUBLE, higherNeighbour,1,_comm, &highRecv);


#endif
}

void UniformPseudoParticleContainer::processTree() {

	int curCellsEdge=1;
	double cellWid[3];

	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<=_maxLevel; curLevel++){
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /=2;
		}

		GatherWellSepLo_Global(cellWid, curCellsEdge, curLevel);

		AllReduceLocalMoments(curCellsEdge, curLevel);

		if(curLevel<_maxLevel) {
			PropagateCellLo_Global(cellWid, curCellsEdge, curLevel);

		}
	}
}
#ifdef ENABLE_MPI
void printTimer(double localTimer, std::string nameOfTimer, MPI_Comm comm){
	double globalTimerMax;
	double globalTimerMin;
	double globalTimerSum;
	int numRanks;
	int myRank;
	MPI_Comm_size(comm,&numRanks);
	MPI_Comm_rank(comm,&myRank);
	MPI_Reduce(&localTimer, &globalTimerMax, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&localTimer, &globalTimerMin, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&localTimer, &globalTimerSum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

	if(myRank == 0){
		std::cout << "\t\t Timer: " << globalTimerMax       		<< "\t\t" <<"s in " << nameOfTimer <<"(Max)" << std::endl;
		std::cout << "\t\t Timer: " << globalTimerMin       		<< "\t\t" <<"s in " << nameOfTimer <<"(Min)" << std::endl;
		std::cout << "\t\t Timer: " << globalTimerSum/numRanks      << "\t\t" <<"s in " << nameOfTimer <<"(Avg)" << std::endl;
	}
}
#endif
void UniformPseudoParticleContainer::printTimers() {
#ifdef ENABLE_MPI

	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_ALL_REDUCE"), "Allreduce", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_GLOBAL"), "CombineMpCellGlobal", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GATHER_WELL_SEP_LO_GLOBAL"), "GatherWellSepLoGlobal", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_GLOBAL"), "PropagateCellLoGlobal", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMBINE_MP_CELL_LOKAL"), "CombineMpCellLokal", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_PROPAGATE_CELL_LO_LOKAL"), "PropagateCellLoLokal", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_COMMUNICATION_HALOS"), "Halo communication", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_HALO_GATHER"), "GatherWellSepLoHalos", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_BUSY_WAITING"), "BusyWaiting", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_FMM_COMPLETE"), "total FMM", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_CALCULATION"), "M2M calculation global", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_INIT"), "M2M Init", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_FINALIZE"), "M2M Finalize", _comm);
	printTimer(global_simulation->timers()->getTime("UNIFORM_PSEUDO_PARTICLE_CONTAINER_GLOBAL_M2M_TRAVERSAL"), "M2M Traversal", _comm);

#endif
}

std::vector<std::vector<MpCell>> &UniformPseudoParticleContainer::getMpCellGlobalTop() {
	return _mpCellGlobalTop;
}

} /* namespace bhfmm */


