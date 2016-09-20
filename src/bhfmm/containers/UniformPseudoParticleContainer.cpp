/*
 * UniformPseudoParticleContainer.cpp
 *
 *  Created on: Feb 5, 2015
 *      Author: tchipevn
 */

#include "UniformPseudoParticleContainer.h"
#include "Simulation.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"
#include "bhfmm/HaloBufferNoOverlap.h"
#include "bhfmm/HaloBufferOverlap.h"
#include <string>
#include <algorithm>
#include <array>

namespace bhfmm {

#ifndef WIGNER
#define WIGNER 0 // 0: original, 1: Wigner rotation acceleration
#endif


#define IsOdd(x) ((x) & 1)
#define ToEven(x) ((x) & ~1)

UniformPseudoParticleContainer::UniformPseudoParticleContainer(
		double domainLength[3], double bBoxMin[3], double bBoxMax[3],
		double LJCellLength[3], unsigned LJSubdivisionFactor, int orderOfExpansions,
		bool periodic) :
		PseudoParticleContainer(orderOfExpansions), _leafContainer(0), _wellSep(1)
		 {

	_periodicBC = periodic;
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
		std::vector<int> neigh = domainDecomp.getNeighbourRanksFullShell();
		int myRank;
//		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
//		std::cout << myRank << " ";
		for(int i = 0; i<26; ++i){
			_neighbours.push_back(neigh[i]);
//			std::cout << neigh[i] << " ";
		}
//		std::cout << "\n";
//		MPI_Barrier(MPI_COMM_WORLD);

	}
	_comm = domainDecomp.getCommunicator();
#endif
#if WIGNER == 1
	//global_log->error() << "not supported yet" << endl;
	exit(-1);
#endif
#ifdef ENABLE_MPI
	_timerProcessCells.set_sync(false);
	_timerAllreduce.set_sync(false);
	_timerAllreduce_me.set_sync(false);
	_timerCombineMpCellGlobal.set_sync(false);
	_timerGatherWellSepLoGlobal.set_sync(false);
	_timerPropagateCellLoGlobal.set_sync(false);
	_timerCombineMpCellLokal.set_sync(false);
	_timerGatherWellSepLoLokal.set_sync(false);
	_timerPropagateCellLoLokal.set_sync(false);
	_timerProcessFarField.set_sync(false);
	_timerCommunicationHalos.set_sync(false);
	_timerHaloGather.set_sync(false);
	_timerBusyWaiting.set_sync(false);
	_timerFMMcomplete.set_sync(false);

#endif
	_leafContainer = new LeafNodesContainer(bBoxMin, bBoxMax, LJCellLength,
			LJSubdivisionFactor, periodic);

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
	assert(_maxLevel == log2(domainLength[1] / cellLength[1]));
	assert(_maxLevel == log2(domainLength[2] / cellLength[2]));
#if defined(ENABLE_MPI)
	int numProcessors;
	MPI_Comm_size(_comm,&numProcessors);
	_globalLevel = ceil(log2(numProcessors)/3.0);
	if(_globalLevel > _maxLevel){
		std::cout << "too many MPI ranks \n";
		exit(-1);
	}
	//numProcessers has to be a power of 2
	assert(pow(2,log2(numProcessors)) == numProcessors);
	//non overlap communication works only for powers of 8
	if(!_overlapComm){
		_numProcessorsPerDim = pow(2,log2(numProcessors)/3);
	}
	for(int i=0; i<3; i++){
		_numProcessorsPerDim[i] = rint(domainLength[i]/(bBoxMax[i] - bBoxMin[i]));
	}
//	std::cout << _numProcessorsPerDim <<"\n\n\n";
	int maxProcsPerDim = std::max(_numProcessorsPerDim[0], std::max(_numProcessorsPerDim[1],_numProcessorsPerDim[2]));
	for(int i=0; i<3; i++){
		_numCellsOnGlobalLevel[i] = rint(maxProcsPerDim * 1.0 /_numProcessorsPerDim[i]);
	}
//	std::cout << _numCellsOnGlobalLevel << _globalLevel <<"\n\n\n";

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
		_mpCellLocal.push_back(std::vector<MpCell>((int) (num_cells_in_level_one_dim[0] + 4) * (num_cells_in_level_one_dim[1] + 4) * (num_cells_in_level_one_dim[2] + 4) , _maxOrd));
		xHaloSize += 2 * num_cells_in_level_one_dim[1] * num_cells_in_level_one_dim[2];
		yHaloSize += 2 * (num_cells_in_level_one_dim[0] + 4) * num_cells_in_level_one_dim[2];
		zHaloSize += 2 * (num_cells_in_level_one_dim[0] + 4) * (num_cells_in_level_one_dim[1] + 4);
		yHaloSizeOverlap += 2 * (num_cells_in_level_one_dim[0]) * num_cells_in_level_one_dim[2];
		zHaloSizeOverlap += 2 * (num_cells_in_level_one_dim[0]) * (num_cells_in_level_one_dim[1]);

		edgeHaloSize += num_cells_in_level_one_dim * 4;
		cornerHaloSize += 8;
		num_cells_in_level_one_dim *= 2;
	}
//	assert(
//			num_cells_in_level
//					== _globalNumCellsPerDim * _globalNumCellsPerDim
//							* _globalNumCellsPerDim);

	// initalize centers and radii
	//num_cells_in_level = 1;
	Vector3<double> current_pos;
	Vector3<double> current_cell_length(domainLength);
	int num_cells_in_level_global = 1;
#if defined(ENABLE_MPI)
	Vector3<double> globalLevelCellLength = Vector3<double>(domainLength)
						* (1.0 / pow(2,_globalLevel));
	_processorPositionGlobalLevel[0] = rint(bBoxMin[0]/ globalLevelCellLength[0]);
	_processorPositionGlobalLevel[1] = rint(bBoxMin[1]/ globalLevelCellLength[1]);
	_processorPositionGlobalLevel[2] = rint(bBoxMin[2]/ globalLevelCellLength[2]);
#endif
//	std::cout <<_processorPositionGlobalLevel << "\n\n";
	//initialization of global top part of tree

	for (int n = 0; n <= _globalLevel; ++n) {
		for (int z = 0; z < num_cells_in_level_global; ++z) {
			for (int y = 0; y < num_cells_in_level_global; ++y) {
				for (int x = 0; x < num_cells_in_level_global; ++x) {
					current_pos[0] = (x + 0.5) * current_cell_length[0];
					current_pos[1] = (y + 0.5) * current_cell_length[1];
					current_pos[2] = (z + 0.5) * current_cell_length[2];
					int cellIndex = ((z * num_cells_in_level_global + y)
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

	//here it is supposed that every processor has exactly one subtree consisting only of one root node
	//ToDo multiple trees from different roots (so far only one)
	//num_cells_in_level = 1;

	//initialization of local tree
//	int num_cells_in_level_one_dim_old = num_cells_in_level_one_dim;
//	num_cells_in_level_one_dim = 2;

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
			//num_cells_in_level *= 8; //ToDo: is it needed anymore?
			current_cell_length = current_cell_length * 0.5;
		}

	_domain = global_simulation->getDomain();
	_globalNumCells = pow(_globalNumCellsPerDim, 3);
#if defined(ENABLE_MPI)
	_globalLevelNumCells = pow(8,_globalLevel);
	_occVector = new int[_globalLevelNumCells];
	std::fill(_occVector, _occVector + _globalLevelNumCells, 0);
#endif
	_coeffVectorLength = 0;
//	_coeffVectorLength = _mpCell[0][0].multipole.get
	for (int j = 0; j <= _maxOrd; j++) {
		for (int k = 0; k <= j; k++) {
			_coeffVectorLength++;
		}
	}
#if defined(ENABLE_MPI)
	xHaloSize *= _coeffVectorLength;
	yHaloSize *= _coeffVectorLength;
	zHaloSize *= _coeffVectorLength;
	yHaloSizeOverlap *= _coeffVectorLength;
	zHaloSizeOverlap *= _coeffVectorLength;
	cornerHaloSize *= _coeffVectorLength;
	edgeHaloSize *= _coeffVectorLength;
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
//		for(int i=0 ; i< 26; i++){
//			std::cout << _neighbours[i] <<" ";
//		}
//		std::cout <<"n";
		_multipoleBufferOverlap = new HaloBufferOverlap<double>(areaHaloSize * 2,edgeHaloSize * 2, cornerHaloSize * 2, _comm, areaNeighbours, edgesNeighbours, cornerNeighbours, 1);
		_multipoleRecBufferOverlap = new HaloBufferOverlap<double>(areaHaloSize * 2,edgeHaloSize * 2, cornerHaloSize * 2, _comm, areaNeighbours, edgesNeighbours, cornerNeighbours, 0);
		_multipoleRecBufferOverlap->communicate();
		//take care that every process has initiated receive in first iteration
		MPI_Barrier(_comm);
	}
#endif

#if defined(ENABLE_MPI)
	int numCells = 0;
	for(int i = _globalLevel; i >=0 ; i--){
		numCells += pow(8,i);
	}
	_coeffVectorLength *= numCells;
#endif

	_coeffVector = new double[_coeffVectorLength * 2];
	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);
	Log::global_log->info() << "UniformPseudoParticleContainer: coeffVectorLength="
			<< _coeffVectorLength << " Size of MPI Buffers is "
			<< (8 * (_coeffVectorLength * 2 + _globalNumCells)
					/ (1024.0 * 1024.0)) << " MB;" << std::endl;
	Log::global_log->info() << "UniformPseudoParticleContainer: globalLevel = "
				<< _globalLevel << " maxLevel= " <<_maxLevel << std::endl;



#ifdef FMM_FFT
	//FFT objects initialization
	FFTSettings::autoSetting(_maxOrd); //auto set FFT acceleration's settings to optimal values
	//FFTSettings::TFMANAGER_VERBOSE = false;
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


}

UniformPseudoParticleContainer::~UniformPseudoParticleContainer() {
	delete _leafContainer;
	delete[] _coeffVector;
#if defined(ENABLE_MPI)
	delete[] _occVector;
	if(!_overlapComm){
		delete _multipoleBuffer;
		delete _multipoleRecBuffer;
	}
	else{
		delete _multipoleBufferOverlap;
		delete _multipoleRecBufferOverlap;
	}
#endif

	
#ifdef FMM_FFT
	if (FFTSettings::issetFFTAcceleration()) {
		delete _FFT_TM;
		delete _FFTAcceleration;
	}
#endif  /* FMM_FFT */
}

void UniformPseudoParticleContainer::build(ParticleContainer* pc) {
	_timerFMMcomplete.start();
	_leafContainer->clearParticles();

	Molecule* tM;
	for(tM = pc->begin(); tM != pc->end(); tM = pc->next()) {
		_leafContainer->addParticle(*tM);
	}
}

void UniformPseudoParticleContainer::upwardPass(P2MCellProcessor* cp) {
	// P2M
	_leafContainer->traverseCells(*cp);

//	if(_maxLevel == _globalLevel){
//		AllReduceMultipoleMoments();
//	}

	_timerCombineMpCellGlobal.start();

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
		if(curLevel >= _globalLevel){
		    Vector3<int> curCellsEdgeLocal;
		    for(int j = 0; j < 3; j++){
		    	curCellsEdgeLocal[j] = (int) (curCellsEdge/_numProcessorsPerDim[j])+4;
		    }
		    const Vector3<int> offset = (_globalLevel == curLevel)? _processorPositionGlobalLevel: Vector3<int>(2);
			CombineMpCell_MPI(cellWid, curCellsEdgeLocal , curLevel, offset);
		}
		else{
			CombineMpCell(cellWid, curCellsEdge, curLevel);
		}
		if(curLevel == _globalLevel + 1){
			communicateHalos();
		}

#else
		CombineMpCell(cellWid, curCellsEdge, curLevel);
#endif

	}
	AllReduceMultipoleMomentsLevelToTop(_globalLevelNumCells,_globalLevel);

	_timerCombineMpCellGlobal.stop();
	//trigger communication by testing with MPI test
	initBusyWaiting();
	busyWaiting();

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
			GatherWellSepLo_FFT_MPI(cellWid, curCellsEdgeLocal, curLevel, 0);
#else
			GatherWellSepLo_MPI(cellWid, curCellsEdgeLocal, curLevel, 0);
#endif
		}
	}
#endif


//	if(_overlapComm){
//		//wait till communication of halos finished
//		_multipoleRecBufferOverlap->wait();
//
//		//set halo values
//		communicateHalosOverlapSetHalos();
//	}
	//wait till allreduce finished
//	MPI_Status allreduceStatus;
//	MPI_Wait(&_allReduceRequest, &allreduceStatus);

	int finishedFlag = 0;
	initBusyWaiting();
	_timerBusyWaiting.start();
	while(finishedFlag != -1){
		finishedFlag = busyWaiting();
#if defined(ENABLE_MPI)

		if(finishedFlag == 1){ //communication of halos finished
			_timerBusyWaiting.stop();
			communicateHalosOverlapSetHalos();
			curCellsEdge=1;
			for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

			for(int curLevel=1; curLevel<=_maxLevel; curLevel++){

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
					GatherWellSepLo_FFT_MPI(cellWid, curCellsEdgeLocal, curLevel, 1);
#else
					GatherWellSepLo_MPI(cellWid, curCellsEdgeLocal, curLevel, 1);
#endif
				}
			}
			//start receiving for next iteration; important for ready send
			_multipoleRecBufferOverlap->communicate();
			_timerBusyWaiting.start();
		}

#endif

		if(finishedFlag == 2){ //communication of global allreduce finished
			_timerBusyWaiting.stop();
			AllReduceMultipoleMomentsSetValues(pow(8,_globalLevel), _globalLevel);
			curCellsEdge=1;
			for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

			for(int curLevel=1; curLevel<=_maxLevel; curLevel++){

				curCellsEdge *=2;
				for(int i=0; i <3; i++){
					cellWid[i] /=2;
				}

#if defined(ENABLE_MPI)
				if(curLevel <= _globalLevel){
	#ifdef FMM_FFT
					GatherWellSepLo_FFT(cellWid, curCellsEdge, curLevel);
	#else
					GatherWellSepLo(cellWid, curCellsEdge, curLevel);
	#endif
				}
#else
	#ifdef FMM_FFT
					GatherWellSepLo_FFT(cellWid, curCellsEdge, curLevel);
	#else
					GatherWellSepLo(cellWid, curCellsEdge, curLevel);
	#endif
					finishedFlag = -1;
#endif
			}
			_timerBusyWaiting.start();
		}
	}
	_timerBusyWaiting.stop();


}

int UniformPseudoParticleContainer::busyWaiting(){
#ifdef ENABLE_MPI
	if(_allReduceProcessed and _halosProcessed and _sendProcessed){
		return -1;
	}
	if(!_sendProcessed){
	    if(_multipoleBufferOverlap->testIfFinished()){
			_sendProcessed = 1;

			return 3;
		}
	}
	if(!_halosProcessed){
		if(_multipoleRecBufferOverlap->testIfFinished()){
			_halosProcessed = 1;

			return 1;
		}
	}
	if(!_allReduceProcessed){
		int flag;
		MPI_Status status;
		MPI_Test(&_allReduceRequest, &flag, &status);
		if(flag){
			_allReduceProcessed = 1;

			return 2;
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


	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<_maxLevel; curLevel++){
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /= 2;
		}

#if defined(ENABLE_MPI) && WIGNER==0
		if(curLevel >= _globalLevel){
			Vector3<int> curCellsEdgeLocal;
			for(int j = 0; j < 3; j++){
				curCellsEdgeLocal[j] = (int) (curCellsEdge/_numProcessorsPerDim[j])+4;
			}
			const Vector3<int> offset = (_globalLevel == curLevel)? _processorPositionGlobalLevel: Vector3<int>(2);
			PropagateCellLo_MPI(cellWid, curCellsEdgeLocal, curLevel,offset);
		}
		else{
			PropagateCellLo(cellWid, curCellsEdge, curLevel);
		}
#else
		PropagateCellLo(cellWid, curCellsEdge, curLevel);
#endif
	}

	// L2P
	_leafContainer->traverseCells(*cp);

	_timerFMMcomplete.stop();
}




void UniformPseudoParticleContainer::CombineMpCell(double */*cellWid*/, int& mpCells, int curLevel){
//#pragma omp parallel firstprivate(curLevel, doHalos, localMpCells)
	{
		int iDir, m1=0, m1x, m1y, m1z, m2=0;
		int m2v[3] = {0, 0, 0};
		int mpCellsN=2*mpCells;
	//#pragma omp for schedule(dynamic,1)
		for (int mloop = 0 ; mloop < mpCells * mpCells * mpCells; mloop++){
			m1x = mloop % mpCells;
			m1y = (mloop / mpCells) % mpCells;
			m1z = mloop / (mpCells * mpCells);
			m1=(m1z*mpCells + m1y)*mpCells + m1x;

			for(iDir=0; iDir<8; iDir++){

				m2v[0]=2*m1x;
				m2v[1]=2*m1y;
				m2v[2]=2*m1z;

				if(IsOdd(iDir  )) m2v[0]=m2v[0]+1;
				if(IsOdd(iDir/2)) m2v[1]=m2v[1]+1;
				if(IsOdd(iDir/4)) m2v[2]=m2v[2]+1;


				m2=(m2v[2]*mpCellsN + m2v[1])*mpCellsN + m2v[0];

				if(_mpCellGlobalTop[curLevel+1][m2].occ==0) continue;

				_mpCellGlobalTop[curLevel][m1].occ +=_mpCellGlobalTop[curLevel+1][m2].occ;

				_mpCellGlobalTop[curLevel][m1].multipole.addMultipoleParticle(_mpCellGlobalTop[curLevel+1][m2].multipole);
			} // iDir closed
		} // mloop closed
	}
}



void UniformPseudoParticleContainer::CombineMpCell_MPI(double* /*cellWid*/, Vector3<int> localMpCells, int curLevel, Vector3<int> offset){
//#pragma omp parallel firstprivate(curLevel, doHalos, localMpCells,offset)
	{
		int iDir, m1=0, m1x, m1y, m1z, m2=0;
		int m2v[3] = {0, 0, 0};
		//take care of halo cells
		Vector3<int> localMpCellsN;
		for(int i = 0; i < 3; i++){
			localMpCellsN[i] = 2 * (localMpCells[i] - 4) + 4;
		}
		Vector3<int> localMpCellsRow;
		std::vector<std::vector<MpCell> > * mpCellCurLevel;
		int curLevelp1 = curLevel -_globalLevel;

		if(curLevel == _globalLevel){
			for(int i = 0; i < 3; i++){
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
		for(int d = 0; d < 3; d++){
			numInnerCells[d] = localMpCells[d] - 4;
		}
		//#pragma omp for schedule(dynamic,1)
		for (int mloop = 0 ; mloop < numInnerCells[0] * numInnerCells[1] * numInnerCells[2]; mloop++){
			m1x = mloop % numInnerCells[0];
			m1y = (mloop / numInnerCells[0]) % numInnerCells[1];
			m1z = mloop / (numInnerCells[0] * numInnerCells[1]);
			m1=((m1z+offset[2])*localMpCellsRow[1] + m1y+offset[1])*localMpCellsRow[0] + m1x+offset[0];

			for(iDir=0; iDir<8; iDir++){

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
}

#define HiLim(t) ToEven(m1v[t])+ 2*_wellSep+1
#define LoLim(t) ToEven(m1v[t])- 2*_wellSep

void UniformPseudoParticleContainer::GatherWellSepLo(double *cellWid, int mpCells, int curLevel){
	_timerGatherWellSepLoGlobal.start();
//#pragma omp parallel firstprivate(curLevel, doHalos, localMpCells)
	{
		int m1v[3];
		int m2v[3];

		int m1, m2, m2x, m2y, m2z;
		// int m1x, m1y, m1z;
		int m22x, m22y, m22z; // for periodic image
		int _size, _rank, loop_min, loop_max;
		int _row_length;

		_row_length=mpCells*mpCells*mpCells;

		DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	//	_rank= domainDecomp.getRank();
	//	_size= domainDecomp.getNumProcs();
	//
	//	loop_min = (int) ((long) (_rank + 0) * (long) (_row_length) / (long) _size);
	//	loop_max = (int) ((long) (_rank + 1) * (long) (_row_length) / (long) _size);

		Vector3<double> periodicShift;
//#pragma omp for schedule(dynamic,1)
		for (m1 = 0; m1 < _row_length; m1++) {

			m1v[0] = m1 % mpCells;
			m1v[1] = (m1 / mpCells) % mpCells;
			m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
			if (_mpCellGlobalTop[curLevel][m1].occ == 0)
				continue;

			for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
				if (_periodicBC == false and (m2z < 0 or m2z >= mpCells)) {
					continue;
				}

				// to get periodic image
				m22z = (mpCells + m2z) % mpCells;
				periodicShift[2] = 0.0;
				if (m2z < 0) 		periodicShift[2] = -mpCells * cellWid[2];
				if (m2z >= mpCells) periodicShift[2] = mpCells * cellWid[2];

				m2v[2] = m2z;
				for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
					if (_periodicBC == false and (m2y < 0 or m2y >= mpCells)) {
						continue;
					}

					// to get periodic image
					m22y = (mpCells + m2y) % mpCells;

					periodicShift[1] = 0.0;
					if (m2y < 0)		periodicShift[1] = -mpCells * cellWid[1];
					if (m2y >= mpCells)	periodicShift[1] = mpCells * cellWid[1];

					m2v[1] = m2y;
					for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
						if (_periodicBC == false and (m2x < 0 or m2x >= mpCells)) {
							continue;
						}

						// to get periodic image
						m22x = (mpCells + m2x) % mpCells;

						periodicShift[0] = 0.0;
						if (m2x < 0)		periodicShift[0] = -mpCells * cellWid[0];
						if (m2x >= mpCells) periodicShift[0] = mpCells * cellWid[0];
						//
						m2v[0] = m2x;

						if (abs(m2v[0] - m1v[0]) <= _wellSep &&
							abs(m2v[1] - m1v[1]) <= _wellSep &&
							abs(m2v[2] - m1v[2]) <= _wellSep)
							continue;
						m2 = (m22z * mpCells + m22y) * mpCells + m22x;
	#ifndef ENABLE_MPI
						if (_mpCellGlobalTop[curLevel][m2].occ == 0)
							continue;
	#endif
						_mpCellGlobalTop[curLevel][m1].local.addMultipoleParticle(
								_mpCellGlobalTop[curLevel][m2].multipole, periodicShift);
					} // m2x closed
				} // m2y closed
			} // m2z closed
		} //m1 closed
	}
	_timerGatherWellSepLoGlobal.stop();

} // GatherWellSepLo closed



void UniformPseudoParticleContainer::GatherWellSepLo_MPI(double* /*cellWid*/, Vector3<int> localMpCells, int curLevel, int doHalos){
	_timerGatherWellSepLoLokal.start();
	if(doHalos){
		_timerHaloGather.start();
	}
//#pragma omp parallel firstprivate(curLevel, doHalos, localMpCells)
	{
		int m1x,m1y,m1z;
		int m2v[3];
		int m1v[3];
		//adjust for local level
		curLevel = curLevel - _globalLevel - 1;
		int m1, m2, m2x, m2y, m2z;

		int _row_length;

		Vector3<double> periodicShift(0.0);
		int zStart = 2;
		int xStart = 2;
		int yStart = 2;
		int xEnd = localMpCells[0] - 2 - xStart;
		int yEnd = localMpCells[1] - 2 - yStart;
		int zEnd = localMpCells[2] - 2 - zStart;


		//#pragma omp for schedule(dynamic,1)
		for (int mloop = 0 ; mloop < xEnd * yEnd * zEnd; mloop++){
			m1x = mloop % xEnd  + xStart;
			m1y = (mloop / xEnd) % yEnd + yStart;
			m1z = mloop / (xEnd * yEnd) + zStart;
			if(doHalos and (m1x >= 4 and m1x >= 4 and m1z >= 4 and m1x < localMpCells[0] - 4 and m1y < localMpCells[1] - 4 and m1z < localMpCells[2] -4)
					//or (m1z < 2 or m1z >= localMpCells[2] - 2 or m1y < 2 or m1y >= localMpCells[1] - 2 or m1x < 2 or m1x >= localMpCells[0] - 2)) //cannot be halo cell
					){ // either inner or halo cell
				continue;
			}

			m1=((m1z)*localMpCells[1] + m1y)*localMpCells[0] + m1x;
			m1v[0] = m1x;
			m1v[1] = m1y;
			m1v[2] = m1z;
			if (_mpCellLocal[curLevel][m1].occ == 0 ){
				continue;
			}
			bool inHaloz, inHaloy, inHalox; //shows if current cell is in halo area in respective coordinate axis
			inHaloz = inHaloy = inHalox = 0;
			for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
				if ((m2z < 2 or m2z >= localMpCells[2] - 2)) { //halo cell
					if(!doHalos){
						continue;
					}
					if(doHalos){
						inHaloz = 1;
					}
				}

				m2v[2] = m2z;
				for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
					if ((m2y < 2 or m2y >= localMpCells[1] - 2)) { //halo cell
						if(!doHalos){
							continue;
						}
						if(doHalos){
							inHaloy = 1;
						}
					}

					m2v[1] = m2y;
					for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
						if ((m2x < 2 or m2x >= localMpCells[0] - 2)) { //halo cell
							if(!doHalos){
								continue;
							}
							if(doHalos){
								inHalox = 1;
							}
						}
						if(doHalos and not(inHalox or inHaloy or inHaloz)){
							continue;
						}
						inHalox = 0;

						m2v[0] = m2x;

						if (abs(m2v[0] - m1x) <= _wellSep &&
							abs(m2v[1] - m1y) <= _wellSep &&
							abs(m2v[2] - m1z) <= _wellSep)
							continue;
						m2 = (m2z * localMpCells[1] + m2y) * localMpCells[0] + m2x;
	//							if(doHalos){
	//								std::cout <<"m2: " <<m2x << m2y << m2z << "\n";
	//							}
						if (_mpCellLocal[curLevel][m2].occ == 0)
							continue;
						_mpCellLocal[curLevel][m1].local.addMultipoleParticle(
								_mpCellLocal[curLevel][m2].multipole, periodicShift);

					} // m2x closed
					inHaloy = 0;
				} // m2y closed
				inHaloz = 0;
			} // m2z closed

		} //mloop closed
	}
	_timerGatherWellSepLoLokal.stop();
	if(doHalos){
		_timerHaloGather.stop();
	}
} // GatherWellSepLo closed


#ifdef FMM_FFT
void UniformPseudoParticleContainer::GatherWellSepLo_FFT(double *cellWid, int mpCells, int& curLevel) {
	if (FFTSettings::USE_VECTORIZATION) {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<true, true, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<true, true, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<true, true, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<true, true, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<true, false, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<true, false, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<true, false, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<true, false, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		}
	} else {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<false, true, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<false, true, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<false, true, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<false, true, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<false, false, true, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<false, false, true, false>(
							cellWid, mpCells, curLevel);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_template<false, false, false, true>(
							cellWid, mpCells, curLevel);
				} else {
					GatherWellSepLo_FFT_template<false, false, false, false>(
							cellWid, mpCells, curLevel);
				}
			}
		}
	}
}
  
void UniformPseudoParticleContainer::GatherWellSepLo_FFT_MPI(double *cellWid, Vector3<int> mpCells, int curLevel, int doHalos) {
	if (FFTSettings::USE_VECTORIZATION) {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_MPI_template<true, true, true, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_MPI_template<true, true, true, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_MPI_template<true, true, false, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_MPI_template<true, true, false, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_MPI_template<true, false, true, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_MPI_template<true, false, true, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_MPI_template<true, false, false, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_MPI_template<true, false, false, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			}
		}
	} else {
		if (FFTSettings::USE_TFMANAGER_UNIFORMGRID) {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_MPI_template<false, true, true, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_MPI_template<false, true, true, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_MPI_template<false, true, false, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_MPI_template<false, true, false, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			}
		} else {
			if (FFTSettings::USE_2WAY_M2L) {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_MPI_template<false, false, true, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_MPI_template<false, false, true, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			} else {
				if (FFTSettings::USE_ORDER_REDUCTION) {
					GatherWellSepLo_FFT_MPI_template<false, false, false, true>(
							cellWid, mpCells, curLevel, doHalos);
				} else {
					GatherWellSepLo_FFT_MPI_template<false, false, false, false>(
							cellWid, mpCells, curLevel, doHalos);
				}
			}
		}
	}
}

template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
void UniformPseudoParticleContainer::GatherWellSepLo_FFT_template(
		double *cellWid, int mpCells, int& curLevel) {
	_timerGatherWellSepLoGlobal.start();
//#pragma omp parallel firstprivate(curLevel, mpCells, cellWid)
	{
		int m1v[3];
		int m2v[3];

		int m1, m2, m2x, m2y, m2z;
		// int m1x, m1y, m1z;
		int m22x, m22y, m22z; // for periodic image
	//	int _size, _rank,
		int loop_min, loop_max;
		int row_length;

		row_length = mpCells * mpCells * mpCells;

	//	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	//	_rank = domainDecomp.getRank();
	//	_size = domainDecomp.getNumProcs();

		loop_min = 0;
		loop_max = row_length;

		Vector3<double> periodicShift;

		//FFT param
		double radius;
		double base_unit = 2.0 / sqrt(3);
		int M2L_order;
		FFTDataContainer* tf;

		//Initialize FFT
//#pragma omp for schedule(dynamic,1)
		for (m1 = loop_min; m1 < loop_max; m1++) {
//			if (_mpCellGlobalTop[curLevel][m1].occ == 0)
//				continue;

			radius = _mpCellGlobalTop[curLevel][m1].local.getRadius();
			FFTAccelerableExpansion& source =
					static_cast<bhfmm::SHMultipoleParticle&>(_mpCellGlobalTop[curLevel][m1].multipole).getExpansion();
			FFTAccelerableExpansion& target =
					static_cast<bhfmm::SHLocalParticle&>(_mpCellGlobalTop[curLevel][m1].local).getExpansion();
			_FFTAcceleration->FFT_initialize_Source(source, radius);
			_FFTAcceleration->FFT_initialize_Target(target);

		}

		//M2L in Fourier space
//#pragma omp for schedule(dynamic,1)
		for (m1 = loop_min; m1 < loop_max; m1++) {

			m1v[0] = m1 % mpCells;
			m1v[1] = (m1 / mpCells) % mpCells;
			m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
			if (_mpCellGlobalTop[curLevel][m1].occ == 0)
				continue;

			for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
				// to get periodic image
				m22z = (mpCells + m2z) % mpCells;
				periodicShift[2] = 0.0;
				if (m2z < 0)
					periodicShift[2] = -mpCells * cellWid[2];
				if (m2z >= mpCells)
					periodicShift[2] = mpCells * cellWid[2];

				m2v[2] = m2z;
				for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
					// to get periodic image
					m22y = (mpCells + m2y) % mpCells;

					periodicShift[1] = 0.0;
					if (m2y < 0)
						periodicShift[1] = -mpCells * cellWid[1];
					if (m2y >= mpCells)
						periodicShift[1] = mpCells * cellWid[1];

					m2v[1] = m2y;
					for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
						// to get periodic image
						m22x = (mpCells + m2x) % mpCells;

						periodicShift[0] = 0.0;
						if (m2x < 0)
							periodicShift[0] = -mpCells * cellWid[0];
						if (m2x >= mpCells)
							periodicShift[0] = mpCells * cellWid[0];
						//
						m2v[0] = m2x;

						if (abs(m2v[0] - m1v[0]) <= _wellSep
								&& abs(m2v[1] - m1v[1]) <= _wellSep
								&& abs(m2v[2] - m1v[2]) <= _wellSep)
							continue;
						m2 = (m22z * mpCells + m22y) * mpCells + m22x;

#ifndef ENABLE_MPI
						if (_mpCellGlobalTop[curLevel][m2].occ == 0){
							continue;
						}
#endif

						if (UseM2L_2way) {
							if (m1 > m2)
								continue;
						}

						tf = _FFT_TM->getTransferFunction(m2x - m1v[0],
								m2y - m1v[1], m2z - m1v[2], base_unit, base_unit,
								base_unit);

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

						if (!UseTFMemoization) {
							delete tf; //free useless memory
						}
						//_mpCell[curLevel][m1].local.addMultipoleParticle(_mpCell[curLevel][m2].multipole, periodicShift);
					} // m2x closed
				} // m2y closed
			} // m2z closed
		} //m1 closed

		//Finalize FFT
//#pragma omp for schedule(dynamic,1)
		for (m1 = loop_min; m1 < loop_max; m1++) {
			if (_mpCellGlobalTop[curLevel][m1].occ == 0)
				continue;

			radius = _mpCellGlobalTop[curLevel][m1].local.getRadius();
			FFTAccelerableExpansion& target =
					static_cast<bhfmm::SHLocalParticle&>(_mpCellGlobalTop[curLevel][m1].local).getExpansion();
			_FFTAcceleration->FFT_finalize_Target(target, radius);
		}
	}
	_timerGatherWellSepLoGlobal.stop();

} // GatherWellSepLo_FFT_template closed

template<bool UseVectorization, bool UseTFMemoization, bool UseM2L_2way, bool UseOrderReduction>
void UniformPseudoParticleContainer::GatherWellSepLo_FFT_MPI_template(
		double *cellWid, Vector3<int> localMpCells, int curLevel, int doHalos) {
	_timerGatherWellSepLoLokal.start();
//#pragma omp parallel firstprivate(curLevel, doHalos, localMpCells)
	{
		//adjust for local level
		curLevel = curLevel - _globalLevel - 1;
		int m1v[3];
		int m2v[3];
		int m1x, m1y, m1z;
		int m1, m2, m2x, m2y, m2z;
		// int m1x, m1y, m1z;
		int m22x, m22y, m22z; // for periodic image
	//	int _size, _rank, loop_min, loop_max;
	//	int row_length;
	//
	//
	//	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	//	_rank = domainDecomp.getRank();
	//	_size = domainDecomp.getNumProcs();


		Vector3<double> periodicShift(0.0);

		//FFT param
		double radius;
		double base_unit = 2.0 / sqrt(3);
		int M2L_order;
		FFTDataContainer* tf;
		int zStart = 2;
		int xStart = 2;
		int yStart = 2;
		int xEnd = localMpCells[0] - 2 - xStart;
		int yEnd = localMpCells[1] - 2 - yStart;
		int zEnd = localMpCells[2] - 2 - zStart;

		if(doHalos){ //in this case also iterate over halo values
			xStart = yStart = zStart = 0;
			xEnd = localMpCells[0];
			yEnd = localMpCells[1];
			zEnd = localMpCells[2];
		}

		//Initialize FFT
//#pragma omp for schedule(dynamic,1)
		for (int mloop = 0 ; mloop < xEnd * yEnd * zEnd; mloop++){
			m1x = mloop % xEnd  + xStart;
			m1y = (mloop / xEnd) % yEnd + yStart;
			m1z = mloop / (xEnd * yEnd) + zStart;

			if(doHalos and not(m1z < 2 or m1z >= localMpCells[2] - 2) and not(m1y < 2 or m1y >= localMpCells[1] - 2) and not(m1x < 2 or m1x >= localMpCells[0] - 2)){
				continue; //only initialize halo values in doHalos run
			}

			m1=((m1z)*localMpCells[1] + m1y)*localMpCells[0] + m1x;
			if (_mpCellLocal[curLevel][m1].occ == 0)
				continue;

			radius = _mpCellLocal[curLevel][m1].local.getRadius();
			FFTAccelerableExpansion& source =
					static_cast<bhfmm::SHMultipoleParticle&>(_mpCellLocal[curLevel][m1].multipole).getExpansion();
			FFTAccelerableExpansion& target =
					static_cast<bhfmm::SHLocalParticle&>(_mpCellLocal[curLevel][m1].local).getExpansion();
			_FFTAcceleration->FFT_initialize_Source(source, radius);
			_FFTAcceleration->FFT_initialize_Target(target);
		}
		int numInnerCells[3];
		for(int d = 0; d < 3; d++){
			numInnerCells[d] = localMpCells[d] - 4;
		}

		//M2L in Fourier space
		//#pragma omp for schedule(dynamic,1)
		for (int mloop = 0 ; mloop < numInnerCells[0] * numInnerCells[1] * numInnerCells[2]; mloop++){
			m1x = mloop % numInnerCells[0]  + 2;
			m1y = (mloop / numInnerCells[0]) % numInnerCells[1] + 2;
			m1z = mloop / (numInnerCells[0] * numInnerCells[1]) + 2;
			if(doHalos and (m1x >= 4 and m1y >= 4 and m1z >= 4 and m1x < localMpCells[0] - 4 and m1y < localMpCells[1] - 4 and m1z < localMpCells[2] -4))
					//or (m1z < 2 or m1z >= localMpCells[2] - 2 or m1y < 2 or m1y >= localMpCells[1] - 2 or m1x < 2 or m1x >= localMpCells[0] - 2))) //cannot be halo cell anymore
			{ // either inner or halo cell
				continue;
			}

			m1=((m1z)*localMpCells[1] + m1y)*localMpCells[0] + m1x;
			m1v[0] = m1x;
			m1v[1] = m1y;
			m1v[2] = m1z;
			if (_mpCellLocal[curLevel][m1].occ == 0 ){
				continue;
			}
			bool inHaloz, inHaloy, inHalox; //shows if current cell is in halo area in respective coordinate axis
			inHaloz = inHaloy = inHalox = 0;
			for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
				if ((m2z < 2 or m2z >= localMpCells[2] - 2)) { //halo cell
					if(!doHalos){
						continue;
					}
					if(doHalos){
						inHaloz = 1;
					}
				}

				m2v[2] = m2z;
				for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
					if ((m2y < 2 or m2y >= localMpCells[1] - 2)) { //halo cell
						if(!doHalos){
							continue;
						}
						if(doHalos){
							inHaloy = 1;
						}
					}

					m2v[1] = m2y;
					for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
						if ((m2x < 2 or m2x >= localMpCells[0] - 2)) { //halo cell
							if(!doHalos){
								continue;
							}
							if(doHalos){
								inHalox = 1;
							}
						}
						if(doHalos and not(inHalox or inHaloy or inHaloz)){
							continue;
						}
						inHalox = 0;

						m2v[0] = m2x;

						if (abs(m2v[0] - m1v[0]) <= _wellSep
								&& abs(m2v[1] - m1v[1]) <= _wellSep
								&& abs(m2v[2] - m1v[2]) <= _wellSep)
							continue;
						m2 = (m2z * localMpCells[1] + m2y) * localMpCells[0] + m2x;

						if (_mpCellLocal[curLevel][m2].occ == 0)
							continue;

						if (UseM2L_2way) {
							if (m1 > m2 && not(doHalos))
								continue;
						}

						tf = _FFT_TM->getTransferFunction(m2v[0] - m1v[0],
								m2v[1] - m1v[1], m2v[2] - m1v[2], base_unit, base_unit,
								base_unit);

						if (UseM2L_2way) {
							FFTAccelerableExpansion& source2 =
									static_cast<bhfmm::SHMultipoleParticle&>(_mpCellLocal[curLevel][m2].multipole).getExpansion();
							FFTAccelerableExpansion& target2 =
									static_cast<bhfmm::SHLocalParticle&>(_mpCellLocal[curLevel][m2].local).getExpansion();
							FFTAccelerableExpansion& source1 =
									static_cast<bhfmm::SHMultipoleParticle&>(_mpCellLocal[curLevel][m1].multipole).getExpansion();
							FFTAccelerableExpansion& target1 =
									static_cast<bhfmm::SHLocalParticle&>(_mpCellLocal[curLevel][m1].local).getExpansion();

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
									static_cast<bhfmm::SHMultipoleParticle&>(_mpCellLocal[curLevel][m2].multipole).getExpansion();
							FFTAccelerableExpansion& target =
									static_cast<bhfmm::SHLocalParticle&>(_mpCellLocal[curLevel][m1].local).getExpansion();

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

						if (!UseTFMemoization) {
							delete tf; //free useless memory
						}
						//_mpCell[curLevel][m1].local.addMultipoleParticle(_mpCell[curLevel][m2].multipole, periodicShift);
					} // m2x closed
					inHaloy = 0;
				} // m2y closed
				inHaloz = 0;
			} // m2z closed

		} //mloop closed

		//Finalize FFT
		if(doHalos){

//#pragma omp for schedule(dynamic,1)
			for (int mloop = 0 ; mloop < numInnerCells[0] * numInnerCells[1] * numInnerCells[2]; mloop++){
				m1x = mloop % numInnerCells[0]  + 2;
				m1y = (mloop / numInnerCells[0]) % numInnerCells[1] + 2;
				m1z = mloop / (numInnerCells[0] * numInnerCells[1]) + 2;

				m1=((m1z)*localMpCells[1] + m1y)*localMpCells[0] + m1x;
				if (_mpCellLocal[curLevel][m1].occ == 0)
					continue;

				radius = _mpCellLocal[curLevel][m1].local.getRadius();
				FFTAccelerableExpansion& target =
						static_cast<bhfmm::SHLocalParticle&>(_mpCellLocal[curLevel][m1].local).getExpansion();
				_FFTAcceleration->FFT_finalize_Target(target, radius);

			}

		}
	}
	_timerGatherWellSepLoLokal.stop();

} // GatherWellSepLo_FFT_MPI_template closed
#endif /* FMM_FFT */


void UniformPseudoParticleContainer::PropagateCellLo(double */*cellWid*/, int mpCells, int curLevel){
	_timerPropagateCellLoGlobal.start();
//#pragma omp parallel firstprivate(curLevel, doHalos, localMpCells)
	{
		int m1v[3];
		int m2v[3];

		int iDir, m1, m1x, m1y, m1z, m2;

		int mpCellsN = 2*mpCells;

	// TODO: parallelization is broken currently, but parallelizing L2L is not all that important
	//	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();


	//	int _rank= domainDecomp.getRank();
	//	int _size= domainDecomp.getNumProcs();

	//	int loop_min = (int) ((long) (_rank + 0) * (long) (mpCells * mpCells * mpCells) / (long) _size);
	//	int loop_max = (int) ((long) (_rank + 1) * (long) (mpCells * mpCells * mpCells) / (long) _size);
		int loop_min = 0;
		int loop_max = mpCells * mpCells * mpCells;
	//#pragma omp for schedule(dynamic,1)
		for (m1 = loop_min; m1 < loop_max; m1++) {

			m1v[0] = m1 % mpCells;
			m1v[1] = (m1 / mpCells) % mpCells;
			m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
			m1x = m1v[0];
			m1y = m1v[1];
			m1z = m1v[2];

			if (_mpCellGlobalTop[curLevel][m1].occ == 0)
				continue;

			for (iDir = 0; iDir < 8; iDir++) {
				m2v[0] = 2 * m1x;
				m2v[1] = 2 * m1y;
				m2v[2] = 2 * m1z;

				if (IsOdd(iDir))     m2v[0] = m2v[0] + 1;
				if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
				if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

				m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];

				_mpCellGlobalTop[curLevel][m1].local.actOnLocalParticle(
						_mpCellGlobalTop[curLevel + 1][m2].local);
			} // iDir
		}
	}
	_timerPropagateCellLoGlobal.stop();
} // PropogateCellLo


void UniformPseudoParticleContainer::PropagateCellLo_MPI(double* /*cellWid*/, Vector3<int> localMpCells, int curLevel, Vector3<int> offset){
	_timerPropagateCellLoLokal.start();
//#pragma omp parallel firstprivate(curLevel, doHalos, localMpCells, offset)
	{
	//	int m1v[3];
		int m2v[3];

		int iDir, m1, m1x, m1y, m1z, m2;

		Vector3<int> localMpCellsN;
		for(int i = 0; i < 3; i++){
			localMpCellsN[i] = 2 * (localMpCells[i] - 4) + 4;
		}
	// TODO: parallelization is broken currently, but parallelizing L2L is not all that important
	//	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
		//correct length in case of globalLevel is reached
		Vector3<int> localMpCellsRow;
		std::vector<std::vector<MpCell> > * mpCellCurLevel;
		int curLevelp1 = curLevel -_globalLevel;

		if(curLevel == _globalLevel){
			for(int i = 0; i < 3; i++){
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
		int numInnerCells[3];
		for(int d = 0; d < 3; d++){
			numInnerCells[d] = localMpCells[d] - 4;
		}
	//#pragma omp for schedule(dynamic,1)
		for (int mloop = 0 ; mloop < numInnerCells[0] * numInnerCells[1] * numInnerCells[2]; mloop++){
			m1x = mloop % numInnerCells[0];
			m1y = (mloop / numInnerCells[0]) % numInnerCells[1];
			m1z = mloop / (numInnerCells[0] * numInnerCells[1]);
			m1=((m1z+offset[2])*localMpCellsRow[1] + m1y+offset[1])*localMpCellsRow[0] + m1x+offset[0];

			if ((*mpCellCurLevel)[curLevel][m1].occ == 0){
				continue;
			}

			for (iDir = 0; iDir < 8; iDir++) {
				//adjust for halo
				m2v[0] = 2 * m1x + 2;
				m2v[1] = 2 * m1y + 2;
				m2v[2] = 2 * m1z + 2;

				if (IsOdd(iDir))     m2v[0] = m2v[0] + 1;
				if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
				if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

				m2 = (m2v[2] * localMpCellsN[1] + m2v[1]) * localMpCellsN[0] + m2v[0];

				(*mpCellCurLevel)[curLevel][m1].local.actOnLocalParticle(
						_mpCellLocal[curLevelp1][m2].local);
			} // iDir
		}

	}

	_timerPropagateCellLoLokal.stop();
} // PropogateCellLo_MPI



void UniformPseudoParticleContainer::processMultipole(ParticleCell& cell){
	int cellIndexV[3];
	std::vector<std::vector<MpCell> > * mpCellMaxLevel;
	int maxLevel;
	if(_maxLevel == _globalLevel){
		mpCellMaxLevel = &_mpCellGlobalTop;
		maxLevel = _maxLevel;
	}
	else{
		mpCellMaxLevel = &_mpCellLocal;
		maxLevel = _maxLevel - _globalLevel - 1;
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
	int cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1])	* _globalNumCellsPerDim) + cellIndexV[0];
#endif

//	assert(cell.isInActiveWindow());

	int currentParticleCount = cell.getMoleculeCount();


	int Occupied = 0;

	// loop over all particles in the cell
	for (int i = 0; i < currentParticleCount; i++) {
		++Occupied;
		Molecule& molecule1 = cell.moleculesAt(i);
		int ni= molecule1.numCharges();

		for(int j=0; j<ni; j++){
			const std::array<double,3> dii = molecule1.charge_d(j);
			const Charge& chargei=static_cast<const Charge&> (molecule1.component()->charge(j));
			double dr[3];

			for(int k=0; k<3; k++){
				dr[k]=molecule1.r(k)+dii[k];
			}	// for k closed

			bhfmm::Vector3<double> site_pos_vec3(dr);
			(*mpCellMaxLevel)[maxLevel][cellIndex].multipole.addSource(site_pos_vec3, chargei.q());

		}// for j closed
	} // current particle closed

	(*mpCellMaxLevel)[maxLevel][cellIndex].occ = Occupied;
}

void UniformPseudoParticleContainer::processFarField(ParticleCell& cell) {
	int cellIndexV[3];
	std::vector<std::vector<MpCell> > * mpCellMaxLevel;

	int maxLevel;
	if(_maxLevel == _globalLevel){
		mpCellMaxLevel = &_mpCellGlobalTop;
		maxLevel = _maxLevel;
	}
	else{
		mpCellMaxLevel = &_mpCellLocal;
		maxLevel = _maxLevel - _globalLevel - 1;
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

	bhfmm::SolidHarmonicsExpansion leLocal(_maxOrd);
	int currentParticleCount = cell.getMoleculeCount();
	double u = 0;
	double uSum = 0.0;
	double f[3] = {0.0, 0.0, 0.0};
	Vector3<double>f_vec3;
	double virialSum=0.0;
	double P_xxSum=0.0;
	double P_yySum=0.0;
	double P_zzSum=0.0;

	// loop over all particles in the cell
	for (int i = 0; i < currentParticleCount; i++) {
		Molecule& molecule1 = cell.moleculesAt(i);
		int ni= molecule1.numCharges();

		for(int j=0; j<ni; j++){
			const std::array<double,3> dii = molecule1.charge_d(j);
			const Charge& chargei=static_cast<const Charge&> (molecule1.component()->charge(j));
			Vector3<double> dr;

			for(int k=0; k<3; k++){
				dr[k]=molecule1.r(k)+dii[k];
			}       // for k closed

			(*mpCellMaxLevel)[maxLevel][cellIndex].local.actOnTarget(dr,chargei.q(),u,f_vec3);
			f[0] = f_vec3[0];
			f[1] = f_vec3[1];
			f[2] = f_vec3[2];

			double virial = 0.0;
			for(int l=0; l<3; l++){
				virial +=-f[l]*dr[l];
			}
			P_xxSum +=0.5*-f[0]*dr[0];
			P_yySum +=0.5*-f[1]*dr[1];
			P_zzSum +=0.5*-f[2]*dr[2];
			molecule1.Fchargeadd(j, f);
			uSum +=0.5*chargei.q()*u;
			virialSum +=0.5*virial;
		}// for j closed
	} // current particle closed

//	_domain->addLocalUpot(uSum);
//	_domain->addLocalVirial(virialSum);
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
	std::fill(_coeffVector, _coeffVector + _coeffVectorLength*2, 0.0);
	std::fill(_occVector, _occVector + _globalLevelNumCells, 0);
#endif
}

void UniformPseudoParticleContainer::AllReduceMultipoleMoments() {
	_timerAllreduce.start();
#ifdef ENABLE_MPI
	int coeffIndex = 0;
	for (int cellIndex = 0; cellIndex < _globalNumCells; cellIndex++) {
		const MpCell & currentCell = _mpCellGlobalTop[_maxLevel][cellIndex];

		// NOTE: coeffIndex modified in following call:
		currentCell.multipole.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

		assert(cellIndex < _globalNumCells);
		_occVector[cellIndex] = currentCell.occ;
	}

	MPI_Allreduce(MPI_IN_PLACE, _coeffVector, _coeffVectorLength*2, MPI_DOUBLE, MPI_SUM, _comm);
	MPI_Allreduce(MPI_IN_PLACE, _occVector, _globalLevelNumCells, MPI_INT, MPI_SUM, _comm);


	coeffIndex = 0;
	for (int cellIndex = 0; cellIndex < _globalLevelNumCells; cellIndex++) {
		MpCell & currentCell = _mpCellGlobalTop[_maxLevel][cellIndex];

		currentCell.occ = _occVector[cellIndex];
		currentCell.multipole.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);

#endif
	_timerAllreduce.stop();
}

void UniformPseudoParticleContainer::AllReduceMultipoleMomentsLevelToTop(int numCellsLevel,int startingLevel) {
	_timerAllreduce.start();
#ifdef ENABLE_MPI

	int coeffIndex = 0;
	int numCellsLevelTemp = numCellsLevel;

	for(int level = startingLevel; level >=0 ; level--){
		for (int cellIndex = 0; cellIndex < numCellsLevelTemp; cellIndex++) {
			const MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];

			// NOTE: coeffIndex modified in following call:
			currentCell.multipole.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

//			assert(cellIndex < numCellsLevelTemp);
//			_occVector[cellIndex] = currentCell.occ;
		}
		numCellsLevelTemp /= 8;
	}

	MPI_Iallreduce(MPI_IN_PLACE, _coeffVector, _coeffVectorLength*2, MPI_DOUBLE, MPI_SUM, _comm, &_allReduceRequest);
	//MPI_Allreduce(MPI_IN_PLACE, _occVector, numCellsLevel, MPI_INT, MPI_SUM, _comm);


#endif
	_timerAllreduce.stop();
}
void UniformPseudoParticleContainer::AllReduceMultipoleMomentsSetValues(int numCellsLevel,int startingLevel) {
#ifdef ENABLE_MPI
	int coeffIndex = 0;
	int numCellsLevelTemp = numCellsLevel;
	numCellsLevelTemp = numCellsLevel;
	for(int level = startingLevel; level >=0 ; level--){
		for (int cellIndex = 0; cellIndex < numCellsLevelTemp; cellIndex++) {

			MpCell & currentCell = _mpCellGlobalTop[level][cellIndex];

//			currentCell.occ = _occVector[cellIndex];
			currentCell.multipole.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

		}
		numCellsLevelTemp /= 8;
	}
	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);
#endif
}
void UniformPseudoParticleContainer::AllReduceLocalMoments(int mpCells, int _curLevel) {
	_timerAllreduce_me.start();

#ifdef ENABLE_MPI

	const int _row_Length=pow(mpCells, 3);
	int coeffIndex = 0;

	coeffIndex = 0;

	for (int cellIndex = 0; cellIndex < _row_Length; cellIndex++) {
		const MpCell & currentCell = _mpCellGlobalTop[_curLevel][cellIndex];

		if(currentCell.occ == 0) continue;
		currentCell.local.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

	}
	MPI_Allreduce(MPI_IN_PLACE, _coeffVector, coeffIndex, MPI_DOUBLE, MPI_SUM, _comm);

	coeffIndex = 0;

	for (int cellIndex = 0; cellIndex < _row_Length; cellIndex++) {
		MpCell & currentCell = _mpCellGlobalTop[_curLevel][cellIndex];

		if(currentCell.occ == 0) continue;
		currentCell.local.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);

#endif
	_timerAllreduce_me.stop();
}


void UniformPseudoParticleContainer::getHaloValues(Vector3<int> localMpCellsBottom,int bottomLevel, double *buffer,
		int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh){
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
					currentCell.multipole.writeValuesToMPIBuffer(buffer, coeffIndex);
				}
			}
		}
		for(int i = 0; i < 3; i++){
			localMpCells[i] = (localMpCells[i] - 4) / 2 + 4;
		}
	}
#endif
}

void UniformPseudoParticleContainer::setHaloValues(Vector3<int> localMpCellsBottom,int bottomLevel, double *bufferRec,
		int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh){
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

					currentCell.multipole.readValuesFromMPIBuffer(bufferRec, coeffIndex);
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
	_timerCommunicationHalos.start();
	if(!_overlapComm){
		communicateHalosNoOverlap();
	}
	else{
		communicateHalosOverlapStart();
//		_multipoleRecBufferOverlap->wait();
//		//_multipoleBufferOverlap->wait();
//		communicateHalosOverlapSetHalos();
	}
	_timerCommunicationHalos.stop();
}
void UniformPseudoParticleContainer::communicateHalosNoOverlap(){

#if defined(ENABLE_MPI)

	_multipoleBuffer->clear();

	Vector3<int> localMpCellsBottom;
	for(int i = 0; i < 3; i++){
		localMpCellsBottom[i] = pow(2,_maxLevel) / _numProcessorsPerDim[i]  + 4;
	}
	//communicate along x axis
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getLeftBuffer(),
			2, 4, 2, -2, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getRightBuffer(),
			-4, -2, 2, -2, 2, -2);

	communicateHalosAlongAxis(_multipoleBuffer->getLeftBuffer(),_multipoleBuffer->getRightBuffer(),_multipoleRecBuffer->getLeftBuffer(),_multipoleRecBuffer->getRightBuffer(),
			_neighbours[0],_neighbours[1],_multipoleBuffer->getXSize()
			);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getLeftBuffer(),
			0, 2, 2, -2, 2, -2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getRightBuffer(),
			-2, 0, 2, -2, 2, -2);

	//communicate along y axis
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getBottomBuffer(),
			0, 0, 2, 4, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getTopBuffer(),
			0, 0, -4, -2, 2, -2);

	communicateHalosAlongAxis(_multipoleBuffer->getBottomBuffer(),_multipoleBuffer->getTopBuffer(),_multipoleRecBuffer->getBottomBuffer(),_multipoleRecBuffer->getTopBuffer(),
			_neighbours[2],_neighbours[3],_multipoleBuffer->getYSize()
			);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getBottomBuffer(),
			0, 0, 0, 2, 2, -2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getTopBuffer(),
			0, 0, -2, 0, 2, -2);

	//communicate along z axis
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getBackBuffer(),
			0, 0, 0, 0, 2, 4);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBuffer->getFrontBuffer(),
			0, 0, 0, 0, -4, -2);

	communicateHalosAlongAxis(_multipoleBuffer->getBackBuffer(),_multipoleBuffer->getFrontBuffer(),_multipoleRecBuffer->getBackBuffer(),_multipoleRecBuffer->getFrontBuffer(),
			_neighbours[4],_neighbours[5],_multipoleBuffer->getZSize()
			);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getBackBuffer(),
			0, 0, 0, 0, 0, 2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBuffer->getFrontBuffer(),
			0, 0, 0, 0, -2, 0);


	MPI_Barrier(_comm);
#endif
}


void UniformPseudoParticleContainer::communicateHalosOverlapStart(){

#if defined(ENABLE_MPI)
	//start receiving
//	_multipoleRecBufferOverlap->startCommunication();
//	_multipoleRecBufferOverlap->communicate();
	//clear buffers
	_multipoleBufferOverlap->clear();

	Vector3<int> localMpCellsBottom;
	for(int i = 0; i < 3; i++){
		localMpCellsBottom[i] = pow(2,_maxLevel) / _numProcessorsPerDim[i]  + 4;
	}
	//fill buffers of the halo areas of the simulation cube

//	int lowDirection[6] = {2,4,2,-2,2,-2};
//	int highDirection[6] = {-4,-2,2,-2,2,-2};
//
//
//	for(int i=6; i < 6; i = i + 2){
//		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[i],
//				lowDirection[(0-i+6)%6], lowDirection[(1-i+6)%6], lowDirection[(2-i+6)%6], lowDirection[(3-i+6)%6], lowDirection[(4-i+6)%6], lowDirection[(5-i+6)%6]);
//
//		getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[i+1],
//				highDirection[(0-i+6)%6], highDirection[(1-i+6)%6], highDirection[(2-i+6)%6], highDirection[(3-i+6)%6], highDirection[(4-i+6)%6], highDirection[(5-i+6)%6]);
//	}

	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[0],
									2, 4, 2, -2, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[1],
									-4, -2, 2, -2, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[2],
									2, -2, 2, 4, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[3],
									2, -2, -4, -2, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[4],
									2, -2, 2, -2, 2, 4);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getAreaBuffers()[5],
									2, -2, 2, -2, -4, -2);
	//fill edges buffers of the halo areas
	//adjacent edges to lower x area

	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[0],
									2, 4, 2, 4, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[2],
									2, 4, -4, -2, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[4],
									2, 4, 2, -2, 2, 4);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[6],
									2, 4, 2, -2, -4, -2);

	//adjacent edges to higher x area
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[1],
									-4, -2, -4, -2, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[3],
									-4, -2, 2, 4, 2, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[5],
									-4, -2, 2, -2, -4, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[7],
									-4, -2, 2, -2, 2, 4);

	//remaining edges
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[8],
									2, -2, 2, 4, 2, 4);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[10],
									2, -2, 2, 4, -4, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[11],
									2, -2, -4, -2, 2, 4);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getEdgeBuffers()[9],
									2, -2, -4, -2, -4, -2);


	//corners

	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[0],
										2, 4, 2, 4, 2, 4);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[1],
										-4, -2, -4, -2, -4, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[2],
										2, 4, 2, 4, -4, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[3],
										-4, -2, -4, -2, 2, 4);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[4],
										2, 4, -4, -2, 2, 4);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[5],
										-4, -2, 2, 4, -4, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[6],
										2, 4, -4, -2, -4, -2);
	getHaloValues(localMpCellsBottom,_maxLevel, _multipoleBufferOverlap->getCornerBuffers()[7],
										-4, -2, 2, 4, 2, 4);

	//start sending
//	_multipoleBufferOverlap->startCommunication();
	_multipoleBufferOverlap->communicate();

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

//	int lowDirection[6] = {0,2,2,-2,2,-2};
//	int highDirection[6] = {-2,0,2,-2,2,-2};
//
//
//	for(int i=0; i<6; i = i + 2){
//		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[i],
//				lowDirection[(0-i+6)%6], lowDirection[(1-i+6)%6], lowDirection[(2-i+6)%6], lowDirection[(3-i+6)%6], lowDirection[(4-i+6)%6], lowDirection[(5-i+6)%6]);
//
//		setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[i+1],
//				highDirection[(0-i+6)%6], highDirection[(1-i+6)%6], highDirection[(2-i+6)%6], highDirection[(3-i+6)%6], highDirection[(4-i+6)%6], highDirection[(5-i+6)%6]);
//	}
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[0],
									0, 2, 2, -2, 2, -2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[1],
									-2, 0, 2, -2, 2, -2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[2],
									2, -2, 0, 2, 2, -2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[3],
									2, -2, -2, 0, 2, -2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[4],
									2, -2, 2, -2, 0, 2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getAreaBuffers()[5],
									2, -2, 2, -2, -2, 0);

	//read edges buffers of the halo areas
	//adjacent edges to lower x area

	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[0],
									0, 2,0,2,2,-2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[2],
									0, 2,-2,0,2,-2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[4],
									0, 2, 2,-2, 0, 2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[6],
									0, 2,2,-2,-2, 0);

	//adjacent edges to higher x area
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[1],
									-2, 0, -2, 0, 2,-2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[3],
									-2, 0, 0, 2, 2,-2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[5],
									-2, 0, 2, -2, -2, 0);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[7],
									-2, 0, 2, -2, 0, 2);

	//remaining edges
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[8],
									2, -2, 0, 2, 0, 2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[10],
									2, -2, 0, 2, -2, 0);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[11],
									2, -2, -2, 0, 0, 2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getEdgeBuffers()[9],
									2, -2, -2, 0, -2, 0);


	//corners

	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[0],
										0, 2, 0, 2, 0, 2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[1],
										-2, 0, -2, 0, -2, 0);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[2],
										0, 2, 0, 2, -2, 0);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[3],
										-2, 0, -2, 0, 0, 2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[4],
										0, 2, -2, 0, 0, 2);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[5],
										-2, 0, 0, 2, -2, 0);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[6],
										0, 2, -2, 0, -2, 0);
	setHaloValues(localMpCellsBottom,_maxLevel, _multipoleRecBufferOverlap->getCornerBuffers()[7],
										-2, 0, 0, 2, 0, 2);

#endif
}


void UniformPseudoParticleContainer::communicateHalosAlongAxis(double * lowerNeighbourBuffer, double * higherNeighbourBuffer,
		double * lowerNeighbourBufferRec, double * higherNeighbourBufferRec,
		int lowerNeighbour, int higherNeighbour, int haloSize
		){
#if defined(ENABLE_MPI)
	MPI_Request low, high;

	MPI_Status lowRecv,highRecv;

	MPI_Isend(lowerNeighbourBuffer, haloSize, MPI_DOUBLE, lowerNeighbour, 1,
			_comm, &low);
	MPI_Isend(higherNeighbourBuffer, haloSize, MPI_DOUBLE, higherNeighbour, 3,
			_comm, &high);

	MPI_Recv(lowerNeighbourBufferRec, haloSize,MPI_DOUBLE, lowerNeighbour,3,_comm, &lowRecv);
	MPI_Recv(higherNeighbourBufferRec, haloSize,MPI_DOUBLE, higherNeighbour,1,_comm, &highRecv);


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

		GatherWellSepLo(cellWid, curCellsEdge, curLevel);

		AllReduceLocalMoments(curCellsEdge, curLevel);

		if(curLevel<_maxLevel) {
			PropagateCellLo(cellWid, curCellsEdge, curLevel);

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
	printTimer(_timerAllreduce.get_etime(),"Allreduce",_comm);

	//std::cout << "\t\t" << _timerAllreduce_me.get_etime()			<< "\t\t" <<"s in Allreduce_me"<<std::endl;
	printTimer(_timerCombineMpCellGlobal.get_etime(),"CombineMpCellGlobal",_comm);

	printTimer(_timerGatherWellSepLoGlobal.get_etime(),"GatherWellSepLoGlobal",_comm);

	printTimer(_timerPropagateCellLoGlobal.get_etime(),"PropagateCellLoGlobal",_comm);

	printTimer(_timerCombineMpCellLokal.get_etime(),"CombineMpCellLokal",_comm);

	printTimer(_timerGatherWellSepLoLokal.get_etime(),"GatherWellSepLoLokal",_comm);

	printTimer(_timerPropagateCellLoLokal.get_etime(),"PropagateCellLoLokal",_comm);

	printTimer(_timerCommunicationHalos.get_etime(),"Halo communication",_comm);

	printTimer(_timerHaloGather.get_etime(),"GatherWellSepLoHalos",_comm);

	printTimer(_timerBusyWaiting.get_etime(),"BusyWaiting",_comm);

	printTimer(_timerFMMcomplete.get_etime(),"total FMM",_comm);
#endif

}


} /* namespace bhfmm */


