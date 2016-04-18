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
#include "bhfmm/utils/RotationParameterLookUp.h"
#include "particleContainer/ParticleContainer.h"

#include <algorithm>

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
		PseudoParticleContainer(orderOfExpansions), _leafContainer(0), _wellSep(1),
		_M2M_Wigner(4, WignerMatrix(orderOfExpansions, true)), _L2L_Wigner(4, WignerMatrix(orderOfExpansions, true)) {

	_periodicBC = periodic;

#if defined(ENABLE_MPI) && defined(NEW_FMM)
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	std::vector<int> neigh = domainDecomp.getNeighbourRanks();
	for(int i = 0; i<6; ++i){
		_neighbours.push_back(neigh[i]);
		//std::cout << neigh[i] << " , ";
	}
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//std::cout << rank<<"\n";

#endif
#if WIGNER == 1 && defined(NEW_FMM)
	//global_log->error() << "not supported yet" << endl;
	exit(-1);
#endif
#ifdef ENABLE_MPI
	_timerProcessCells.set_sync(false);
	_timerAllreduce.set_sync(false);
	_timerAllreduce_me.set_sync(false);
	_timerCombineMpCell.set_sync(false);
	_timerGatherWellSepLo.set_sync(false);
	_timerPropagateCellLo.set_sync(false);
	_timerProcessFarField.set_sync(false);
#endif
	//ToDo adjust everything to local sizes
	_leafContainer = new LeafNodesContainer(bBoxMin, bBoxMax, LJCellLength,
			LJSubdivisionFactor, periodic);

	double cellLength[3];

	for (int i = 0; i < 3; i++) {
		cellLength[i] = _leafContainer->getCellLength()[i];
		_cellLength[i] = cellLength[i];
	}
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	_bBoxMin[0] = bBoxMin[0];
	_bBoxMin[1] = bBoxMin[1];
	_bBoxMin[2] = bBoxMin[2];
#endif

	_globalNumCellsPerDim = domainLength[0] / cellLength[0];
	_maxLevel = log2(_globalNumCellsPerDim);
	assert(_maxLevel == log2(domainLength[1] / cellLength[1]));
	assert(_maxLevel == log2(domainLength[2] / cellLength[2]));
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	int numProcessors;
	MPI_Comm_size(MPI_COMM_WORLD,&numProcessors);
	_globalLevel = log2(numProcessors)/3;
	if(_globalLevel > _maxLevel){
		std::cout << "too many MPI ranks \n";
		exit(-1);
	}
	//numProcessers has to be a power of 8
	assert(log2(numProcessors) == _globalLevel * 3);
	_numProcessorsPerDim = pow(2,log2(numProcessors)/3);
#endif


	//allocate Multipole and Local particles
	int num_cells_in_level = 1;
	_mpCell.reserve(_maxLevel + 1);
	int num_cells_in_level_one_dim;
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	num_cells_in_level_one_dim = 1;

	for (int n = 0; n <= _globalLevel; n++) {
		_mpCell.push_back(std::vector<MpCell>(num_cells_in_level, _maxOrd));
		num_cells_in_level *= 8;
	}
	//num_cells_in_level = 8;
	num_cells_in_level_one_dim = 2;
	_xHaloSize = 0;
	_yHaloSize = 0;
	_zHaloSize = 0;
	for (int n = _globalLevel + 1; n <= _maxLevel; n++) {
		_mpCell.push_back(std::vector<MpCell>((int) pow(num_cells_in_level_one_dim + 4 , 3), _maxOrd));
		_xHaloSize += 2 * num_cells_in_level_one_dim * num_cells_in_level_one_dim;
		_yHaloSize += 2 * (num_cells_in_level_one_dim + 4) * num_cells_in_level_one_dim;
		_zHaloSize += 2 * (num_cells_in_level_one_dim + 4) * (num_cells_in_level_one_dim + 4);

		num_cells_in_level_one_dim *= 2;
	}
#else
	for (int n = 0; n <= _maxLevel; n++) {
		_mpCell.push_back(std::vector<MpCell>(num_cells_in_level, _maxOrd));
		num_cells_in_level *= 8;
	}
#endif
//	assert(
//			num_cells_in_level
//					== _globalNumCellsPerDim * _globalNumCellsPerDim
//							* _globalNumCellsPerDim);

	// initalize centers and radii
	num_cells_in_level = 1;
	num_cells_in_level_one_dim = 1;
	Vector3<double> current_pos;
	Vector3<double> current_cell_length(domainLength);
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	Vector3<double> globalLevelCellLength = Vector3<double>(domainLength)
						* (1.0 / pow(2,_globalLevel));
	_processorPositionGlobalLevel[0] = bBoxMin[0]/ globalLevelCellLength[0];
	_processorPositionGlobalLevel[1] = bBoxMin[1]/ globalLevelCellLength[1];
	_processorPositionGlobalLevel[2] = bBoxMin[2]/ globalLevelCellLength[2];

	for (int n = 0; n <= _globalLevel; ++n) {
		for (int z = 0; z < num_cells_in_level_one_dim; ++z) {
			for (int y = 0; y < num_cells_in_level_one_dim; ++y) {
				for (int x = 0; x < num_cells_in_level_one_dim; ++x) {
					current_pos[0] = (x + 0.5) * current_cell_length[0];
					current_pos[1] = (y + 0.5) * current_cell_length[1];
					current_pos[2] = (z + 0.5) * current_cell_length[2];
					int cellIndex = ((z * num_cells_in_level_one_dim + y)
							* num_cells_in_level_one_dim) + x;
					_mpCell[n][cellIndex].multipole.setCenter(current_pos);
					_mpCell[n][cellIndex].multipole.setRadius(
							current_cell_length.L2Norm() * 0.5);

					_mpCell[n][cellIndex].local.setCenter(current_pos);
					_mpCell[n][cellIndex].local.setRadius(
							current_cell_length.L2Norm() * 0.5);
				}
			}
		}
		num_cells_in_level_one_dim *= 2;
		num_cells_in_level *= 8; //ToDo: is it needed anymore?
		current_cell_length = Vector3<double>(domainLength)
				* (1.0 / num_cells_in_level_one_dim);
	}

//	ToDo check if coincides with domain decomposition localX,localY,localZ, _globalLevel
//	int localX,localY,localZ;
	//here it is supposed that every processor has exactly one subtree consisting only of one root node
	//ToDo multiple trees from different roots (so far only one)
	//num_cells_in_level = 1;
	int num_cells_in_level_one_dim_old = num_cells_in_level_one_dim;
	num_cells_in_level_one_dim = 2;
//	int xPosition, yPosition, zPosition;
//	int cellsPerGlobalDimension = pow(2,_globalLevel);
//	int myRank;
//	MPI_Comm_rank(MPI_IN_PLACE,&myRank);
//	localX = myRank % cellsPerGlobalDimension;
//	localY = (myRank / cellsPerGlobalDimension) % cellsPerGlobalDimension;
//	localZ = (myRank / (cellsPerGlobalDimension * cellsPerGlobalDimension)) % cellsPerGlobalDimension;
	for (int n = _globalLevel+1; n <= _maxLevel; ++n) {
			for (int z = -2; z < num_cells_in_level_one_dim+2; ++z) {
				for (int y = -2; y < num_cells_in_level_one_dim+2; ++y) {
					for (int x = -2; x < num_cells_in_level_one_dim+2; ++x) {
//						zPosition = z + localZ;
//						yPosition = y + localY;
//						xPosition = x + localX;
						current_pos[0] = (x + 0.5) * current_cell_length[0] + bBoxMin[0];
						current_pos[1] = (y + 0.5) * current_cell_length[1] + bBoxMin[1];
						current_pos[2] = (z + 0.5) * current_cell_length[2] + bBoxMin[2];
						int cellIndex = (((z + 2) * (num_cells_in_level_one_dim + 4) + y + 2)
								* (num_cells_in_level_one_dim + 4)) + x + 2;
						_mpCell[n][cellIndex].multipole.setCenter(current_pos);
						_mpCell[n][cellIndex].multipole.setRadius(
								current_cell_length.L2Norm() * 0.5);

						_mpCell[n][cellIndex].local.setCenter(current_pos);
						_mpCell[n][cellIndex].local.setRadius(
								current_cell_length.L2Norm() * 0.5);
					}
				}
			}
			num_cells_in_level_one_dim *= 2;
			//num_cells_in_level *= 8; //ToDo: is it needed anymore?
			current_cell_length = Vector3<double>(domainLength)
					* (1.0 / (num_cells_in_level_one_dim * _numProcessorsPerDim));
		}
#else
	for (int n = 0; n <= _maxLevel; ++n) {
		for (int z = 0; z < num_cells_in_level_one_dim; ++z) {
			for (int y = 0; y < num_cells_in_level_one_dim; ++y) {
				for (int x = 0; x < num_cells_in_level_one_dim; ++x) {
					current_pos[0] = (x + 0.5) * current_cell_length[0];
					current_pos[1] = (y + 0.5) * current_cell_length[1];
					current_pos[2] = (z + 0.5) * current_cell_length[2];
					int cellIndex = ((z * num_cells_in_level_one_dim + y)
							* num_cells_in_level_one_dim) + x;
					_mpCell[n][cellIndex].multipole.setCenter(current_pos);
					_mpCell[n][cellIndex].multipole.setRadius(
							current_cell_length.L2Norm() * 0.5);

					_mpCell[n][cellIndex].local.setCenter(current_pos);
					_mpCell[n][cellIndex].local.setRadius(
							current_cell_length.L2Norm() * 0.5);
				}
			}
		}
		num_cells_in_level_one_dim *= 2;
		num_cells_in_level *= 8; //ToDo: is it needed anymore?
		current_cell_length = Vector3<double>(domainLength)
				* (1.0 / num_cells_in_level_one_dim);
	}
#endif
	_domain = global_simulation->getDomain();
	_globalNumCells = pow(_globalNumCellsPerDim, 3);
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	_globalLevelNumCells = pow(8,_globalLevel);
	_occVector = new int[_globalLevelNumCells];
	std::fill(_occVector, _occVector + _globalLevelNumCells, 0);
#else
	_occVector = new int[_globalNumCells];
	std::fill(_occVector, _occVector + _globalNumCells, 0);
#endif
	_coeffVectorLength = 0;
//	_coeffVectorLength = _mpCell[0][0].multipole.get
	for (int j = 0; j <= _maxOrd; j++) {
		for (int k = 0; k <= j; k++) {
			_coeffVectorLength++;
		}
	}
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	_xHaloOccSize = _xHaloSize;
	_yHaloOccSize = _yHaloSize;
	_zHaloOccSize = _zHaloSize;

	_xHaloSize *= _coeffVectorLength;
	_yHaloSize *= _coeffVectorLength;
	_zHaloSize *= _coeffVectorLength;

	_leftOccBuffer = new int[_xHaloOccSize];
	_rightOccBuffer = new int[_xHaloOccSize];

	_bottomOccBuffer = new int[_yHaloOccSize];
	_topOccBuffer = new int[_yHaloOccSize];

	_backOccBuffer = new int[_zHaloOccSize];
	_frontOccBuffer = new int[_zHaloOccSize];

	_leftOccBufferRec = new int[_xHaloOccSize];
	_rightOccBufferRec = new int[_xHaloOccSize];

	_bottomOccBufferRec = new int[_yHaloOccSize];
	_topOccBufferRec = new int[_yHaloOccSize];

	_backOccBufferRec = new int[_zHaloOccSize];
	_frontOccBufferRec = new int[_zHaloOccSize];

	_leftBuffer = new double[_xHaloSize * 2];
	_rightBuffer = new double[_xHaloSize * 2];

	_bottomBuffer = new double[_yHaloSize * 2];
	_topBuffer = new double[_yHaloSize * 2];

	_backBuffer = new double[_zHaloSize * 2];
	_frontBuffer = new double[_zHaloSize * 2];

	_leftBufferRec = new double[_xHaloSize * 2];
	_rightBufferRec = new double[_xHaloSize * 2];

	_bottomBufferRec = new double[_yHaloSize * 2];
	_topBufferRec = new double[_yHaloSize * 2];

	_backBufferRec = new double[_zHaloSize * 2];
	_frontBufferRec = new double[_zHaloSize * 2];
#endif

#if defined(ENABLE_MPI) && defined(NEW_FMM)
	_coeffVectorLength *= _globalLevelNumCells;
#else
	_coeffVectorLength *= _globalNumCells;
#endif

	_coeffVector = new double[_coeffVectorLength * 2];
	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);
	Log::global_log->info() << "UniformPseudoParticleContainer: coeffVectorLength="
			<< _coeffVectorLength << " Size of MPI Buffers is "
			<< (8 * (_coeffVectorLength * 2 + _globalNumCells)
					/ (1024.0 * 1024.0)) << " MB;" << std::endl;

	/* Initialize sin/cos(phi) lookup */
	double phi[4];
	phi[0] = -3*M_PI/4;
	phi[1] = -M_PI/4;
	phi[2] = 3*M_PI/4;
	phi[3] = M_PI/4;

	_CosSin = new double[(_maxOrd+1)*4];
	// We only use negative phi -> cos(-x)=cos(x), sin(-x)=-sin(x)
	for (int m = 0; m <= _maxOrd; ++m) {
		CosSin_ptr(0)[2*m] 		= cos(m*phi[0]);
		CosSin_ptr(0)[2*m + 1] 	= sin(m*phi[0]);
		CosSin_ptr(1)[2*m] 		= cos(m*phi[1]);
		CosSin_ptr(1)[2*m + 1] 	= sin(m*phi[1]);
	}

	/* Initialize M2M-Wigner matrices */
	double theta = atan(sqrt(2.0));
	_M2M_Wigner[0].evaluate(M_PI-theta);
	_M2M_Wigner[1].evaluate(theta);
	_M2M_Wigner[2].evaluate(-(M_PI-theta));
	_M2M_Wigner[3].evaluate(-theta);
	/* Initialize Lookup table */
	RotationParameterLookUp::tab = new RotationParameterLookUp(_maxOrd);
	RotationParameterLookUp::tab->initFromDirEval();

	// pre-multiply prefactors
	for (unsigned i = 0; i<_M2M_Wigner.size(); ++i) {
		WignerMatrix &W = _M2M_Wigner[i];
		for (int l = 0; l <= _maxOrd; ++l) {
			for (int m = 0; m <= l; ++m) {
				for (int k = -l; k<=l; ++k) {
					const double factor = RotationParameterLookUp::tab->acc_c(l,m,k);
					W.acc(l,m,k) *= factor ;
				}
			}
		}

	}

	_L2L_Wigner[0].evaluate(M_PI-theta);
	_L2L_Wigner[1].evaluate(theta);
	_L2L_Wigner[2].evaluate(-(M_PI-theta));
	_L2L_Wigner[3].evaluate(-theta);

	// pre-multiply prefactors
	for (unsigned i = 0; i<_L2L_Wigner.size(); ++i) {
		WignerMatrix &W = _L2L_Wigner[i];
		for (int l = 0; l <= _maxOrd; ++l) {
			for (int m = 0; m <= l; ++m) {
				for (int k = -l; k<=l; ++k) {
					const double factor = RotationParameterLookUp::tab->acc_c(l,k,m);
					W.acc(l,m,k) *= factor ;
				}
			}
		}

	}

	/* Initialize M2L-Wigner matrices */
	processTreeInitM2LWigner();

	delete RotationParameterLookUp::tab;


}

UniformPseudoParticleContainer::~UniformPseudoParticleContainer() {
	delete _leafContainer;
	delete[] _coeffVector;
	delete[] _occVector;
	delete _CosSin;
#if defined(NEW_FMM)
	delete[] _leftBuffer;
	delete[] _rightBuffer;
	delete[] _topBuffer;
	delete[] _bottomBuffer;
	delete[] _frontBuffer;
	delete[] _backBuffer;
	delete[] _leftBufferRec;
	delete[] _rightBufferRec;
	delete[] _topBufferRec;
	delete[] _bottomBufferRec;
	delete[] _frontBufferRec;
	delete[] _backBufferRec;
	delete[] _leftOccBuffer;
	delete[] _rightOccBuffer;
	delete[] _topOccBuffer;
	delete[] _bottomOccBuffer;
	delete[] _frontOccBuffer;
	delete[] _backOccBuffer;
	delete[] _leftOccBufferRec;
	delete[] _rightOccBufferRec;
	delete[] _topOccBufferRec;
	delete[] _bottomOccBufferRec;
	delete[] _frontOccBufferRec;
	delete[] _backOccBufferRec;
#endif

}

void UniformPseudoParticleContainer::build(ParticleContainer* pc) {
	_leafContainer->clearParticles();

	Molecule* tM;
	for(tM = pc->begin(); tM != pc->end(); tM = pc->next()) {
		_leafContainer->addParticle(*tM);
	}
}

void UniformPseudoParticleContainer::upwardPass(P2MCellProcessor* cp) {
	// P2M
	_leafContainer->traverseCellPairs(*cp);
#if not defined(NEW_FMM)
	// M2M
	AllReduceMultipoleMoments();
#endif
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	if(_maxLevel == _globalLevel){
		AllReduceMultipoleMoments();
	}
#endif

	_timerCombineMpCell.start();

	int curCellsEdge=_globalNumCellsPerDim;
	double cellWid[3];

	for(int i=0; i <3; i++)	cellWid[i]=_cellLength[i];

	// when considering periodic boundary conditions, there is actually work up to level 1!
	for(int curLevel=_maxLevel-1; curLevel>=1; curLevel--){
		//ToDo adjust curCellsEdge for local mpCells

		curCellsEdge /=2;
		for(int i=0; i <3; i++)	cellWid[i] *=2;
#if defined(ENABLE_MPI) && defined(NEW_FMM)
		if(curLevel >= _globalLevel){
		    int curCellsEdgeLocal = (int) (curCellsEdge/_numProcessorsPerDim)+4;
		    const Vector3<int> offset = (_globalLevel == curLevel)? _processorPositionGlobalLevel: Vector3<int>(2);
		    //std::cout << "Level " <<curLevel << " offset" << offset << "\n";
			CombineMpCell_MPI(cellWid, curCellsEdgeLocal , curLevel, offset);
		}
		else{
			CombineMpCell(cellWid, curCellsEdge, curLevel);
		}
		if(curLevel == _globalLevel){
			AllReduceMultipoleMomentsLevel(_globalLevelNumCells,curLevel);
		}
#elif WIGNER==0
		CombineMpCell(cellWid, curCellsEdge, curLevel);
#else
		CombineMpCell_Wigner(cellWid, curCellsEdge, curLevel);
#endif
	}
	_timerCombineMpCell.stop();
}

void UniformPseudoParticleContainer::horizontalPass(
		VectorizedChargeP2PCellProcessor* cp) {
	// P2P
	_leafContainer->traverseCellPairs(*cp);

	// M2L
	int curCellsEdge=1;
	double cellWid[3];

	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<=_maxLevel; curLevel++){
		//ToDo adjust curCellsEdge for local mpCells

		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /=2;
		}

#if defined(ENABLE_MPI) && defined(NEW_FMM)
		if(curLevel > _globalLevel){
		    int curCellsEdgeLocal = (int) (curCellsEdge/_numProcessorsPerDim)+4;
			GatherWellSepLo_MPI(cellWid, curCellsEdgeLocal, curLevel);
		}
		else{
			GatherWellSepLo(cellWid, curCellsEdge, curLevel);
		}
#elif WIGNER==0
		GatherWellSepLo(cellWid, curCellsEdge, curLevel);
#else
		GatherWellSepLo_Wigner(cellWid, curCellsEdge, curLevel);
#endif
#if defined(ENABLE_MPI) && defined(NEW_FMM)
		if(curLevel <= _globalLevel){
			AllReduceLocalMoments(curCellsEdge, curLevel);
		}
#else
		AllReduceLocalMoments(curCellsEdge, curLevel);
#endif
	}
}

void UniformPseudoParticleContainer::downwardPass(L2PCellProcessor* cp) {
	// L2L
	int curCellsEdge=1;
	double cellWid[3];


	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<_maxLevel; curLevel++){
		//ToDo adjust curCellsEdge for local mpCells
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /= 2;
		}

#if defined(ENABLE_MPI) && defined(NEW_FMM)
		if(curLevel >= _globalLevel){
		    int curCellsEdgeLocal = (int) (curCellsEdge/_numProcessorsPerDim)+4;
		    const Vector3<int> offset = (_globalLevel == curLevel)? _processorPositionGlobalLevel: Vector3<int>(2);
		    //std::cout << curLevel << " " << offset << " " << curCellsEdgeLocal << "\n";
			PropagateCellLo_MPI(cellWid, curCellsEdgeLocal, curLevel,offset);
		}
		else{
			PropagateCellLo(cellWid, curCellsEdge, curLevel);
		}
#elif WIGNER==0
		PropagateCellLo(cellWid, curCellsEdge, curLevel);
#else
		PropagateCellLo_Wigner(cellWid, curCellsEdge, curLevel);
#endif
	}

	// L2P
	_leafContainer->traverseCellPairs(*cp);
}



void UniformPseudoParticleContainer::CombineMpCell(double *cellWid, int& mpCells, int& curLevel){
	int iDir, m1=0, m1x, m1y, m1z, m2=0;
	int m2v[3] = {0, 0, 0};
	int mpCellsN=2*mpCells;


	for(m1z=0; m1z<mpCells; m1z++){
		for(m1y=0; m1y<mpCells; m1y++){
			for(m1x=0; m1x<mpCells; m1x++){
				m1=(m1z*mpCells + m1y)*mpCells + m1x;

				for(iDir=0; iDir<8; iDir++){

					m2v[0]=2*m1x;
					m2v[1]=2*m1y;
					m2v[2]=2*m1z;

					if(IsOdd(iDir  )) m2v[0]=m2v[0]+1;
					if(IsOdd(iDir/2)) m2v[1]=m2v[1]+1;
					if(IsOdd(iDir/4)) m2v[2]=m2v[2]+1;


					m2=(m2v[2]*mpCellsN + m2v[1])*mpCellsN + m2v[0];

					if(_mpCell[curLevel+1][m2].occ==0) continue;

					_mpCell[curLevel][m1].occ +=_mpCell[curLevel+1][m2].occ;

					_mpCell[curLevel][m1].multipole.addMultipoleParticle(_mpCell[curLevel+1][m2].multipole);
				} // iDir closed
			}// m1x closed
		}// m1y closed
	} // m1z closed
}

void UniformPseudoParticleContainer::CombineMpCell_Wigner(double *cellWid, int& mpCells, int& curLevel){
	int iDir, m1=0, m1x, m1y, m1z, m2=0;
	int m2v[3] = {0, 0, 0};
	int mpCellsN=2*mpCells;

	const double magnitude = sqrt(3.0)/4.0 *cellWid[0];

	for(m1z=0; m1z<mpCells; m1z++){
		for(m1y=0; m1y<mpCells; m1y++){
			for(m1x=0; m1x<mpCells; m1x++){
				m1=(m1z*mpCells + m1y)*mpCells + m1x;

				for(iDir=0; iDir<8; iDir++){

					m2v[0]=2*m1x;
					m2v[1]=2*m1y;
					m2v[2]=2*m1z;

					if(IsOdd(iDir  )) m2v[0]=m2v[0]+1;
					if(IsOdd(iDir/2)) m2v[1]=m2v[1]+1;
					if(IsOdd(iDir/4)) m2v[2]=m2v[2]+1;


					m2=(m2v[2]*mpCellsN + m2v[1])*mpCellsN + m2v[0];

					if(_mpCell[curLevel+1][m2].occ==0) continue;

					_mpCell[curLevel][m1].occ +=_mpCell[curLevel+1][m2].occ;

					/* get Wigner matrix */
					const int idxPhi = iDir%2;
					const int negatePhi = 1 - ((iDir%4)/2)*2;
					const int idxWig = iDir/4;

					_mpCell[curLevel][m1].multipole.addMultipoleParticle_Wigner(_mpCell[curLevel+1][m2].multipole, M2M_Wigner(idxWig),
							M2M_Wigner(idxWig+2), CosSin_ptr(idxPhi), negatePhi, magnitude);
				} // iDir closed
			}// m1x closed
		}// m1y closed
	} // m1z closed
}

void UniformPseudoParticleContainer::CombineMpCell_MPI(double *cellWid, int& localMpCells, int& curLevel, Vector3<int> offset){
	int iDir, m1=0, m1x, m1y, m1z, m2=0;
	int m2v[3] = {0, 0, 0};
	//take care of halo cells
	int localMpCellsN=2*(localMpCells-4)+4;
	int localMpCellsRow;
	if(curLevel == _globalLevel){
		localMpCellsRow = (localMpCells-4) * _numProcessorsPerDim;
	}
	else{
		localMpCellsRow = localMpCells;
	}

	for(m1z=0; m1z<localMpCells-4; m1z++){
		for(m1y=0; m1y<localMpCells-4; m1y++){
			for(m1x=0; m1x<localMpCells-4; m1x++){
				m1=((m1z+offset[2])*localMpCellsRow + m1y+offset[1])*localMpCellsRow + m1x+offset[0];

				for(iDir=0; iDir<8; iDir++){

					m2v[0]=2*m1x+2;
					m2v[1]=2*m1y+2;
					m2v[2]=2*m1z+2;

					if(IsOdd(iDir  )) m2v[0]=m2v[0]+1;
					if(IsOdd(iDir/2)) m2v[1]=m2v[1]+1;
					if(IsOdd(iDir/4)) m2v[2]=m2v[2]+1;


					m2=(m2v[2]*localMpCellsN + m2v[1])*localMpCellsN + m2v[0];

					if(_mpCell[curLevel+1][m2].occ==0) continue;

					_mpCell[curLevel][m1].occ +=_mpCell[curLevel+1][m2].occ;

					_mpCell[curLevel][m1].multipole.addMultipoleParticle(_mpCell[curLevel+1][m2].multipole);
				} // iDir closed
			}// m1x closed
		}// m1y closed
	} // m1z closed
}

#define HiLim(t) ToEven(m1v[t])+ 2*_wellSep+1
#define LoLim(t) ToEven(m1v[t])- 2*_wellSep

void UniformPseudoParticleContainer::GatherWellSepLo(double *cellWid, int mpCells, int& curLevel){
	_timerGatherWellSepLo.start();

	int m1v[3];
	int m2v[3];

	int m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
	int m22x, m22y, m22z; // for periodic image
	int _size, _rank, loop_min, loop_max;
	int _row_length;

	_row_length=mpCells*mpCells*mpCells;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_rank= domainDecomp.getRank();
	_size= domainDecomp.getNumProcs();

	loop_min = (int) ((long) (_rank + 0) * (long) (_row_length) / (long) _size);
	loop_max = (int) ((long) (_rank + 1) * (long) (_row_length) / (long) _size);

	Vector3<double> periodicShift;

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		if (_mpCell[curLevel][m1].occ == 0)
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

					if (_mpCell[curLevel][m2].occ == 0)
						continue;

					_mpCell[curLevel][m1].local.addMultipoleParticle(
							_mpCell[curLevel][m2].multipole, periodicShift);
				} // m2x closed
			} // m2y closed
		} // m2z closed
	} //m1 closed

	_timerGatherWellSepLo.stop();

} // GatherWellSepLo closed

void UniformPseudoParticleContainer::GatherWellSepLo_Wigner(double *cellWid, int mpCells, int& curLevel){
	_timerGatherWellSepLo.start();

	int m1v[3];
	int m2v[3];

	int m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
	int m22x, m22y, m22z; // for periodic image
	int _size, _rank, loop_min, loop_max;
	int _row_length;

	_row_length=mpCells*mpCells*mpCells;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_rank= domainDecomp.getRank();
	_size= domainDecomp.getNumProcs();

	loop_min = (int) ((long) (_rank + 0) * (long) (_row_length) / (long) _size);
	loop_max = (int) ((long) (_rank + 1) * (long) (_row_length) / (long) _size);

	Vector3<double> periodicShift;

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		if (_mpCell[curLevel][m1].occ == 0)
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

					if (_mpCell[curLevel][m2].occ == 0)
						continue;

					_mpCell[curLevel][m1].local.addMultipoleParticle_Wigner(
							_mpCell[curLevel][m2].multipole, periodicShift, cellWid, _M2L_Wigner);
				} // m2x closed
			} // m2y closed
		} // m2z closed
	} //m1 closed

	_timerGatherWellSepLo.stop();

}

void UniformPseudoParticleContainer::GatherWellSepLo_MPI(double *cellWid, int localMpCells, int& curLevel){
	_timerGatherWellSepLo.start();

	int m1x,m1y,m1z;
	int m2v[3];
	int m1v[3];


	int m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
//	int m22x, m22y, m22z; // for periodic image
//	int _size, _rank, loop_min, loop_max;
	int _row_length;

//	_row_length=localMpCells*localMpCells*localMpCells;

//	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
//	_rank= domainDecomp.getRank();
//	_size= domainDecomp.getNumProcs();

//	loop_min = localMpCells*localMpCells+2;
//	loop_max = _row_length-(localMpCells*localMpCells+2);

	Vector3<double> periodicShift(0.0);

	for (m1z = 2; m1z < localMpCells-2; m1z++) {
		for (m1y = 2; m1y < localMpCells-2; m1y++) {
			for (m1x = 2; m1x < localMpCells-2; m1x++) {
				m1=((m1z)*localMpCells + m1y)*localMpCells + m1x;
				m1v[0] = m1x - 2;
				m1v[1] = m1y - 2;
				m1v[2] = m1z - 2;
//		m1v[0] = m1 % localMpCells;
//		m1v[1] = (m1 / localMpCells) % localMpCells;
//		m1v[2] = (m1 / (localMpCells * localMpCells)) % localMpCells;
				if (_mpCell[curLevel][m1].occ == 0 ){
					continue;
				}

				for (m2z = LoLim(2) + 2; m2z <= HiLim(2) + 2; m2z++) {
					if (m2z < 0 or m2z >= localMpCells) {
						std::cout << "Error \n";
						exit(-1);
					}
//
//			// to get periodic image
//			m22z = (mpCells + m2z) % mpCells;
//			periodicShift[2] = 0.0;
//			if (m2z < 0) 		periodicShift[2] = -mpCells * cellWid[2];
//			if (m2z >= mpCells) periodicShift[2] = mpCells * cellWid[2];

					m2v[2] = m2z;
					for (m2y = LoLim(1) + 2; m2y <= HiLim(1) + 2; m2y++) {
						if (m2y < 0 or m2y >= localMpCells) {
							std::cout << "Error \n";
							exit(-1);
						}
//
//				// to get periodic image
//				m22y = (mpCells + m2y) % mpCells;
//
//				periodicShift[1] = 0.0;
//				if (m2y < 0)		periodicShift[1] = -mpCells * cellWid[1];
//				if (m2y >= mpCells)	periodicShift[1] = mpCells * cellWid[1];

						m2v[1] = m2y;
						for (m2x = LoLim(0) + 2; m2x <= HiLim(0) + 2; m2x++) {
							if (m2x < 0 or m2x >= localMpCells) {
								std::cout << "Error \n";
								exit(-1);
							}
//
//					// to get periodic image
//					m22x = (mpCells + m2x) % mpCells;
//
//					periodicShift[0] = 0.0;
//					if (m2x < 0)		periodicShift[0] = -mpCells * cellWid[0];
//					if (m2x >= mpCells) periodicShift[0] = mpCells * cellWid[0];
//					//
							m2v[0] = m2x;

							if (abs(m2v[0] - m1x) <= _wellSep &&
								abs(m2v[1] - m1y) <= _wellSep &&
								abs(m2v[2] - m1z) <= _wellSep)
								continue;
							m2 = (m2z * localMpCells + m2y) * localMpCells + m2x;

							if (_mpCell[curLevel][m2].occ == 0)
								continue;

							_mpCell[curLevel][m1].local.addMultipoleParticle(
									_mpCell[curLevel][m2].multipole, periodicShift);
						} // m2x closed
					} // m2y closed
				} // m2z closed
			} //m1x closed
		} //m1y closed
	} //m1z closed

	_timerGatherWellSepLo.stop();

} // GatherWellSepLo closed

void UniformPseudoParticleContainer::PropagateCellLo(double *cellWid, int mpCells, int& curLevel){
	_timerPropagateCellLo.start();
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

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		m1x = m1v[0];
		m1y = m1v[1];
		m1z = m1v[2];

		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		for (iDir = 0; iDir < 8; iDir++) {
			m2v[0] = 2 * m1x;
			m2v[1] = 2 * m1y;
			m2v[2] = 2 * m1z;

			if (IsOdd(iDir))     m2v[0] = m2v[0] + 1;
			if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
			if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

			m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];

			_mpCell[curLevel][m1].local.actOnLocalParticle(
					_mpCell[curLevel + 1][m2].local);
		} // iDir
	}
	_timerPropagateCellLo.stop();
} // PropogateCellLo

void UniformPseudoParticleContainer::PropagateCellLo_Wigner(double *cellWid, int mpCells, int& curLevel){
	_timerPropagateCellLo.start();
	int m1v[3];
	int m2v[3];

	int iDir, m1, m1x, m1y, m1z, m2;

	int mpCellsN = 2*mpCells;
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	int _rank= domainDecomp.getRank();
	int _size= domainDecomp.getNumProcs();

	int loop_min = (int) ((long) (_rank + 0) * (long) (mpCells * mpCells * mpCells) / (long) _size);
	int loop_max = (int) ((long) (_rank + 1) * (long) (mpCells * mpCells * mpCells) / (long) _size);

	const double magnitude = sqrt(3.0)/4.0 *cellWid[0];

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
		m1x = m1v[0];
		m1y = m1v[1];
		m1z = m1v[2];

		if (_mpCell[curLevel][m1].occ == 0)
			continue;

		for (iDir = 0; iDir < 8; iDir++) {
			m2v[0] = 2 * m1x;
			m2v[1] = 2 * m1y;
			m2v[2] = 2 * m1z;

			if (IsOdd(iDir))     m2v[0] = m2v[0] + 1;
			if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
			if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

			m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];

			/* get Wigner matrix */
			const int idxPhi = iDir%2;
			const int negatePhi = 1 - ((iDir%4)/2)*2;
			const int idxWig = iDir/4;

			_mpCell[curLevel][m1].local.actOnLocalParticle_Wigner(
					_mpCell[curLevel + 1][m2].local, L2L_Wigner(idxWig),
							L2L_Wigner(idxWig+2), CosSin_ptr(idxPhi), negatePhi, magnitude);

		} // iDir
	}
	_timerPropagateCellLo.stop();
} // PropogateCellLo

void UniformPseudoParticleContainer::PropagateCellLo_MPI(double *cellWid, int localMpCells, int& curLevel, Vector3<int> offset){
	_timerPropagateCellLo.start();
//	int m1v[3];
	int m2v[3];

	int iDir, m1, m1x, m1y, m1z, m2;

	int mpCellsN = 2*(localMpCells-4) + 4;

// TODO: parallelization is broken currently, but parallelizing L2L is not all that important
//	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	//correct length in case of globalLevel is reached
	int localMpCellsRow;
	if(curLevel == _globalLevel){
		localMpCellsRow = (localMpCells-4) * _numProcessorsPerDim;
	}
	else{
		localMpCellsRow = localMpCells;
	}
//	int _rank= domainDecomp.getRank();
//	int _size= domainDecomp.getNumProcs();

//	int loop_min = (int) ((long) (_rank + 0) * (long) (mpCells * mpCells * mpCells) / (long) _size);
//	int loop_max = (int) ((long) (_rank + 1) * (long) (mpCells * mpCells * mpCells) / (long) _size);
//	int loop_min = localMpCells * localMpCells + 2;
//	int loop_max = localMpCells * localMpCells * localMpCells - (localMpCells * localMpCells + 2);
	for (m1z = 0; m1z < localMpCells-4; m1z++) {
		for (m1y = 0; m1y < localMpCells-4; m1y++) {
			for (m1x = 0; m1x < localMpCells-4; m1x++) {
				m1=((m1z+offset[2])*localMpCellsRow + m1y+offset[1])*localMpCellsRow + m1x+offset[0];
				if (_mpCell[curLevel][m1].occ == 0)
					continue;

				for (iDir = 0; iDir < 8; iDir++) {
					//adjust for halo
					m2v[0] = 2 * m1x + 2;
					m2v[1] = 2 * m1y + 2;
					m2v[2] = 2 * m1z + 2;

					if (IsOdd(iDir))     m2v[0] = m2v[0] + 1;
					if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
					if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;

					m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];

					_mpCell[curLevel][m1].local.actOnLocalParticle(
							_mpCell[curLevel + 1][m2].local);
				} // iDir
			}
		}
	}
//		m1v[0] = m1 % localMpCells;
//		m1v[1] = (m1 / localMpCells) % localMpCells;
//		m1v[2] = (m1 / (localMpCells * localMpCells)) % localMpCells;
//		//adjust for halo
//		m1x = m1v[0] - 2;
//		m1y = m1v[1] - 2;
//		m1z = m1v[2] - 2;

//		if (_mpCell[curLevel][m1].occ == 0)
//			continue;
//
//		for (iDir = 0; iDir < 8; iDir++) {
//			//adjust for halo
//			m2v[0] = 2 * m1x + 2;
//			m2v[1] = 2 * m1y + 2;
//			m2v[2] = 2 * m1z + 2;
//
//			if (IsOdd(iDir))     m2v[0] = m2v[0] + 1;
//			if (IsOdd(iDir / 2)) m2v[1] = m2v[1] + 1;
//			if (IsOdd(iDir / 4)) m2v[2] = m2v[2] + 1;
//
//			m2 = (m2v[2] * mpCellsN + m2v[1]) * mpCellsN + m2v[0];
//
//			_mpCell[curLevel][m1].local.actOnLocalParticle(
//					_mpCell[curLevel + 1][m2].local);
//		} // iDir

	_timerPropagateCellLo.stop();
} // PropogateCellLo_MPI



void UniformPseudoParticleContainer::processMultipole(ParticleCell& cell){
	int cellIndexV[3];
	//std::cout << "processMultipole: " << cell.getBoxMin(0)<< " "  << cell.getBoxMin(1)<< " "  << cell.getBoxMin(2) << " " <<_bBoxMin << "\n";
	for (int i = 0; i < 3; i++) {
#if defined(ENABLE_MPI) && defined(NEW_FMM)
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
	//std::cout << "ci " << cellIndexV[0] <<" , " <<cellIndexV[1] << " , " <<cellIndexV[2] << " \n";
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	int cellIndex;
	if (_maxLevel == _globalLevel){
		cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1]) * _globalNumCellsPerDim) + cellIndexV[0];
	}
	else{
		const int numCellsOnLocalMaxLevel = (_globalNumCellsPerDim / _numProcessorsPerDim + 4);
		cellIndex = ((numCellsOnLocalMaxLevel * cellIndexV[2] + cellIndexV[1])	* numCellsOnLocalMaxLevel) + cellIndexV[0];
	}
#else
	int cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1])	* _globalNumCellsPerDim) + cellIndexV[0];
#endif

//	assert(cell.isInActiveWindow());

	std::vector<Molecule*>& currentCellParticles = cell.getParticlePointers();
	int currentParticleCount = currentCellParticles.size();


	int Occupied = 0;

	// loop over all particles in the cell
	for (int i = 0; i < currentParticleCount; i++) {
		++Occupied;
		Molecule& molecule1 = *currentCellParticles[i];
		int ni= molecule1.numCharges();

		for(int j=0; j<ni; j++){
			const double* dii = molecule1.charge_d(j);
			const Charge& chargei=static_cast<const Charge&> (molecule1.component()->charge(j));
			double dr[3];

			for(int k=0; k<3; k++){
				dr[k]=molecule1.r(k)+dii[k];
			}	// for k closed

			bhfmm::Vector3<double> site_pos_vec3(dr);
			_mpCell[_maxLevel][cellIndex].multipole.addSource(site_pos_vec3, chargei.q());

		}// for j closed
	} // current particle closed

	_mpCell[_maxLevel][cellIndex].occ = Occupied;
}

void UniformPseudoParticleContainer::processFarField(ParticleCell& cell) {
	int cellIndexV[3];
	for (int i = 0; i < 3; i++) {
	#if defined(ENABLE_MPI) && defined(NEW_FMM)
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

#if defined(ENABLE_MPI) && defined(NEW_FMM)
	int cellIndex;
	if (_maxLevel == _globalLevel){
		cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1]) * _globalNumCellsPerDim) + cellIndexV[0];
	}
	else{
		const int numCellsOnLocalMaxLevel = (_globalNumCellsPerDim / _numProcessorsPerDim + 4);
		cellIndex = ((numCellsOnLocalMaxLevel * cellIndexV[2] + cellIndexV[1])	* numCellsOnLocalMaxLevel) + cellIndexV[0];
	}
#else
	int cellIndex = ((_globalNumCellsPerDim * cellIndexV[2] + cellIndexV[1]) * _globalNumCellsPerDim) + cellIndexV[0];
#endif

	bhfmm::SolidHarmonicsExpansion leLocal(_maxOrd);
	std::vector<Molecule*>& currentCellParticles = cell.getParticlePointers();
	int currentParticleCount = currentCellParticles.size();
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
		Molecule& molecule1 = *currentCellParticles[i];
		int ni= molecule1.numCharges();

		for(int j=0; j<ni; j++){
			const double* dii = molecule1.charge_d(j);
			const Charge& chargei=static_cast<const Charge&> (molecule1.component()->charge(j));
			Vector3<double> dr;

			for(int k=0; k<3; k++){
				dr[k]=molecule1.r(k)+dii[k];
			}       // for k closed

			_mpCell[_maxLevel][cellIndex].local.actOnTarget(dr,chargei.q(),u,f_vec3);
			f[0] = f_vec3[0];
			f[1] = f_vec3[1];
			f[2] = f_vec3[2];
			//if(f[0] > 1 || f[1] > 1 || f[2] >1)
				//std::cout << f[0] << " , " << f[1] << " , " << f[2] << "\n";
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

#if defined(ENABLE_MPI) && defined(NEW_FMM)
	for (int n = _maxLevel; n > _globalLevel; n--) {
		int localMpCells = pow(2, n) / _numProcessorsPerDim + 4;

		for (int m1z = 0; m1z < localMpCells; m1z++) {
			for (int m1y = 0; m1y < localMpCells; m1y++) {
				for (int m1x = 0; m1x < localMpCells; m1x++) {
					int cellIndexNew = (m1z * localMpCells + m1y) * localMpCells + m1x;

					_mpCell[n][cellIndexNew].occ = 0;

					_mpCell[n][cellIndexNew].multipole.clear();
					_mpCell[n][cellIndexNew].local.clear();
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

					_mpCell[n][cellIndexNew].occ = 0;

					_mpCell[n][cellIndexNew].multipole.clear();
					_mpCell[n][cellIndexNew].local.clear();
				}
			}
		}
	}
#else
	for (int n = _maxLevel; n >= 1; n--) {
		int mpCells = pow(2, n);

		for (int m1z = 0; m1z < mpCells; m1z++) {
			for (int m1y = 0; m1y < mpCells; m1y++) {
				for (int m1x = 0; m1x < mpCells; m1x++) {
					int cellIndexNew = (m1z * mpCells + m1y) * mpCells + m1x;

					_mpCell[n][cellIndexNew].occ = 0;

					_mpCell[n][cellIndexNew].multipole.clear();
					_mpCell[n][cellIndexNew].local.clear();
				}
			}
		}
	}
#endif
	// clear the MPI buffers
#ifdef ENABLE_MPI
	std::fill(_coeffVector, _coeffVector + _coeffVectorLength*2, 0.0);
#ifdef NEW_FMM
	std::fill(_occVector, _occVector + _globalLevelNumCells, 0);

#else
	std::fill(_occVector, _occVector + _globalNumCells, 0);
#endif
#endif
}

void UniformPseudoParticleContainer::AllReduceMultipoleMoments() {
	_timerAllreduce.start();
#ifdef ENABLE_MPI
	//std::cout << "error \n";
	int coeffIndex = 0;
	for (int cellIndex = 0; cellIndex < _globalNumCells; cellIndex++) {
		const MpCell & currentCell = _mpCell[_maxLevel][cellIndex];

		// NOTE: coeffIndex modified in following call:
		currentCell.multipole.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

		assert(cellIndex < _globalNumCells);
		_occVector[cellIndex] = currentCell.occ;
	}

	MPI_Allreduce(MPI_IN_PLACE, _coeffVector, _coeffVectorLength*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifdef NEW_FMM
	MPI_Allreduce(MPI_IN_PLACE, _occVector, _globalLevelNumCells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	MPI_Allreduce(MPI_IN_PLACE, _occVector, _globalNumCells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#endif

	coeffIndex = 0;
#ifdef NEW_FMM
	for (int cellIndex = 0; cellIndex < _globalLevelNumCells; cellIndex++) {
#else
	for (int cellIndex = 0; cellIndex < _globalNumCells; cellIndex++) {
#endif
		MpCell & currentCell = _mpCell[_maxLevel][cellIndex];

		currentCell.occ = _occVector[cellIndex];
		currentCell.multipole.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);

#endif
	_timerAllreduce.stop();
}

void UniformPseudoParticleContainer::AllReduceMultipoleMomentsLevel(int numCellsLevel,int curLevel) {
	_timerAllreduce.start();
#ifdef ENABLE_MPI

	int coeffIndex = 0;
	//std::cout << "mcl index: " <<numCellsLevel<< " "<< _globalLevel <<"\n";

	for (int cellIndex = 0; cellIndex < numCellsLevel; cellIndex++) {
		const MpCell & currentCell = _mpCell[curLevel][cellIndex];

		// NOTE: coeffIndex modified in following call:
		currentCell.multipole.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

		assert(cellIndex < numCellsLevel);
		_occVector[cellIndex] = currentCell.occ;
	}

	MPI_Allreduce(MPI_IN_PLACE, _coeffVector, _coeffVectorLength*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, _occVector, numCellsLevel, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	coeffIndex = 0;
	for (int cellIndex = 0; cellIndex < numCellsLevel; cellIndex++) {

		MpCell & currentCell = _mpCell[curLevel][cellIndex];

		currentCell.occ = _occVector[cellIndex];
		currentCell.multipole.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);

#endif
	_timerAllreduce.stop();
}

void UniformPseudoParticleContainer::AllReduceLocalMoments(int mpCells, int _curLevel) {
	_timerAllreduce_me.start();

#ifdef ENABLE_MPI

	const int _row_Length=pow(mpCells, 3);
	int coeffIndex = 0;

	coeffIndex = 0;

	for (int cellIndex = 0; cellIndex < _row_Length; cellIndex++) {
		const MpCell & currentCell = _mpCell[_curLevel][cellIndex];

		if(currentCell.occ == 0) continue;
		currentCell.local.writeValuesToMPIBuffer(_coeffVector, coeffIndex);

	}
	//std::cout << "co index: " <<coeffIndex << " " << _row_Length <<"\n";
	MPI_Allreduce(MPI_IN_PLACE, _coeffVector, coeffIndex, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	coeffIndex = 0;

	for (int cellIndex = 0; cellIndex < _row_Length; cellIndex++) {
		MpCell & currentCell = _mpCell[_curLevel][cellIndex];

		if(currentCell.occ == 0) continue;
		currentCell.local.readValuesFromMPIBuffer(_coeffVector, coeffIndex);

	}

	std::fill(_coeffVector, _coeffVector + _coeffVectorLength * 2, 0.0);

#endif
	_timerAllreduce_me.stop();
}

void UniformPseudoParticleContainer::getXHaloValues(int localMpCellsBottom,int bottomLevel){
#if defined(ENABLE_MPI) && defined(NEW_FMM)

	int coeffIndex = 0;
	int localMpCells = localMpCellsBottom;
	int cellIndex;
	int coeffOccIndex = 0;
	for(int level=bottomLevel; level>_globalLevel;level--){

		//left Border
		for (int z = 2; z < localMpCells-2; z++) {
			for (int y = 2; y < localMpCells-2; y++) {
				for (int x = 2; x < 4; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					const MpCell & currentCell = _mpCell[level][cellIndex];

					_leftOccBuffer[coeffOccIndex] = currentCell.occ;
					coeffOccIndex++;

					//if(currentCell.occ == 0) continue;
					currentCell.multipole.writeValuesToMPIBuffer(_leftBuffer, coeffIndex);
					//std::cout << "Indices: "<< coeffIndex << " " << coeffOccIndex <<"\n";
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
	coeffIndex = 0;
	coeffOccIndex = 0;
	localMpCells = localMpCellsBottom;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//right Border
		for (int z = 2; z < localMpCells-2; z++) {
			for (int y = 2; y < localMpCells-2; y++) {
				for (int x = localMpCells-4; x < localMpCells-2; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					const MpCell & currentCell = _mpCell[level][cellIndex];

					_rightOccBuffer[coeffOccIndex] = currentCell.occ;
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.writeValuesToMPIBuffer(_rightBuffer, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
#endif
}

void UniformPseudoParticleContainer::getYHaloValues(int localMpCellsBottom,int bottomLevel){
#if defined(ENABLE_MPI) && defined(NEW_FMM)

	int coeffIndex = 0;
	int localMpCells = localMpCellsBottom;
	int cellIndex;
	int coeffOccIndex = 0;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//bottom Border
		for (int z = 2; z < localMpCells-2; z++) {
			for (int y = 2; y < 4; y++) {
				for (int x = 0; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					const MpCell & currentCell = _mpCell[level][cellIndex];

					_bottomOccBuffer[coeffOccIndex] = currentCell.occ;
					coeffOccIndex++;

					//if(currentCell.occ == 0) continue;
					currentCell.multipole.writeValuesToMPIBuffer(_bottomBuffer, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
	coeffIndex = 0;
	coeffOccIndex = 0;
	localMpCells = localMpCellsBottom;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//top Border
		for (int z = 2; z < localMpCells-2; z++) {
			for (int y = localMpCells-4; y < localMpCells-2; y++) {
				for (int x = 0; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					const MpCell & currentCell = _mpCell[level][cellIndex];

					_topOccBuffer[coeffOccIndex] = currentCell.occ;
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.writeValuesToMPIBuffer(_topBuffer, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
#endif
}

void UniformPseudoParticleContainer::getZHaloValues(int localMpCellsBottom,int bottomLevel){
#if defined(ENABLE_MPI) && defined(NEW_FMM)

	int coeffIndex = 0;
	int localMpCells = localMpCellsBottom;
	int cellIndex;
	int coeffOccIndex = 0;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//back Border
		for (int z = 2; z < 4; z++) {
			for (int y = 0; y < localMpCells; y++) {
				for (int x = 0; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					const MpCell & currentCell = _mpCell[level][cellIndex];

					_backOccBuffer[coeffOccIndex] = currentCell.occ;
					coeffOccIndex++;

					//if(currentCell.occ == 0) continue;
					currentCell.multipole.writeValuesToMPIBuffer(_backBuffer, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
	coeffIndex = 0;
	coeffOccIndex = 0;
	localMpCells = localMpCellsBottom;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//front Border
		for (int z = localMpCells-4; z < localMpCells-2; z++) {
			for (int y = 0; y < localMpCells; y++) {
				for (int x = 0; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					const MpCell & currentCell = _mpCell[level][cellIndex];

					_frontOccBuffer[coeffOccIndex] = currentCell.occ;
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.writeValuesToMPIBuffer(_frontBuffer, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
#endif
}

void UniformPseudoParticleContainer::setXHaloValues(int localMpCellsBottom,int bottomLevel){
#if defined(ENABLE_MPI) && defined(NEW_FMM)

	int coeffIndex = 0;
	int coeffOccIndex = 0;
	int localMpCells = localMpCellsBottom;
	int cellIndex;
	//std::cout << bottomLevel<< " " << _globalLevel <<";\n";
	//std::cout << localMpCells<< " MPcells   " << _numProcessorsPerDim <<";\n";

	for(int level=bottomLevel; level>_globalLevel;level--){
		//left Border
		for (int z = 2; z < localMpCells-2; z++) {
			for (int y = 2; y < localMpCells-2; y++) {
				for (int x = 0; x < 2; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					MpCell & currentCell = _mpCell[level][cellIndex];
					//std::cout << "before " <<currentCell.multipole.getCenter();

					currentCell.occ = _leftOccBufferRec[coeffOccIndex];
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.readValuesFromMPIBuffer(_leftBufferRec, coeffIndex);
					//std::cout << "after " <<currentCell.multipole.getCenter() << "\n";

				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
	coeffIndex = 0;
	coeffOccIndex = 0;
	localMpCells = localMpCellsBottom;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//right Border
		for (int z = 2; z < localMpCells-2; z++) {
			for (int y = 2; y < localMpCells-2; y++) {
				for (int x = localMpCells-2; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					MpCell & currentCell = _mpCell[level][cellIndex];

					currentCell.occ = _rightOccBufferRec[coeffOccIndex];
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.readValuesFromMPIBuffer(_rightBufferRec, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
#endif
}

void UniformPseudoParticleContainer::setYHaloValues(int localMpCellsBottom,int bottomLevel){
#if defined(ENABLE_MPI) && defined(NEW_FMM)

	int coeffIndex = 0;
	int localMpCells = localMpCellsBottom;
	int cellIndex;
	int coeffOccIndex = 0;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//bottom Border
		for (int z = 2; z < localMpCells-2; z++) {
			for (int y = 0; y < 2; y++) {
				for (int x = 0; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					MpCell & currentCell = _mpCell[level][cellIndex];

					currentCell.occ = _bottomOccBufferRec[coeffOccIndex];
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.readValuesFromMPIBuffer(_bottomBufferRec, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
	coeffIndex = 0;
	coeffOccIndex = 0;
	localMpCells = localMpCellsBottom;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//top Border
		for (int z = 2; z < localMpCells-2; z++) {
			for (int y = localMpCells-2; y < localMpCells; y++) {
				for (int x = 0; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					MpCell & currentCell = _mpCell[level][cellIndex];

					currentCell.occ = _topOccBufferRec[coeffOccIndex];
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.readValuesFromMPIBuffer(_topBufferRec, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
#endif
}

void UniformPseudoParticleContainer::setZHaloValues(int localMpCellsBottom,int bottomLevel){
#if defined(ENABLE_MPI) && defined(NEW_FMM)

	int coeffIndex = 0;
	int coeffOccIndex = 0;

	int localMpCells = localMpCellsBottom;
	int cellIndex;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//back Border
		for (int z = 0; z < 2; z++) {
			for (int y = 0; y < localMpCells; y++) {
				for (int x = 0; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					MpCell & currentCell = _mpCell[level][cellIndex];

					currentCell.occ = _backOccBufferRec[coeffOccIndex];
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.readValuesFromMPIBuffer(_backBufferRec, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
	coeffIndex = 0;
	coeffOccIndex = 0;

	localMpCells = localMpCellsBottom;
	for(int level=bottomLevel; level>_globalLevel;level--){
		//front Border
		for (int z = localMpCells-2; z < localMpCells; z++) {
			for (int y = 0; y < localMpCells; y++) {
				for (int x = 0; x < localMpCells; x++) {
					cellIndex = (z * localMpCells + y) * localMpCells + x;
					MpCell & currentCell = _mpCell[level][cellIndex];

					currentCell.occ = _frontOccBufferRec[coeffOccIndex];
					coeffOccIndex++;
					//if(currentCell.occ == 0) continue;
					currentCell.multipole.readValuesFromMPIBuffer(_frontBufferRec, coeffIndex);
				}
			}
		}
		localMpCells = (localMpCells - 4) / 2 + 4;
	}
#endif
}

void UniformPseudoParticleContainer::communicateHalos(){

#if defined(ENABLE_MPI) && defined(NEW_FMM)
	std::fill(_leftBuffer, _leftBuffer + _xHaloSize * 2, 0.0);
	std::fill(_rightBuffer, _rightBuffer + _xHaloSize * 2, 0.0);
	std::fill(_frontBuffer, _frontBuffer + _zHaloSize * 2, 0.0);
	std::fill(_backBuffer, _backBuffer + _zHaloSize * 2, 0.0);
	std::fill(_topBuffer, _topBuffer + _yHaloSize * 2, 0.0);
	std::fill(_bottomBuffer, _bottomBuffer + _yHaloSize * 2, 0.0);

	//std::cout << "Sizes: "<< _xHaloSize<< " " << _yHaloSize <<" " << _zHaloSize << "\n" ;
	int localMpCellsBottom = pow(2,_maxLevel) / _numProcessorsPerDim  + 4;
	getXHaloValues(localMpCellsBottom,_maxLevel);
	communicateHalosX();
	setXHaloValues(localMpCellsBottom,_maxLevel);

	getYHaloValues(localMpCellsBottom,_maxLevel);
	communicateHalosY();
	setYHaloValues(localMpCellsBottom,_maxLevel);

	getZHaloValues(localMpCellsBottom,_maxLevel);
	communicateHalosZ();
	setZHaloValues(localMpCellsBottom,_maxLevel);
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}
void UniformPseudoParticleContainer::communicateHalosX(){
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	MPI_Request low, high;
	MPI_Request lowOcc, highOcc;

	MPI_Status lowRecv,highRecv;
	MPI_Status lowOccRecv,highOccRecv;

	MPI_Isend(_leftBuffer, _xHaloSize * 2, MPI_DOUBLE, _neighbours[0], 1,
			MPI_COMM_WORLD, &low);
	MPI_Isend(_rightBuffer, _xHaloSize * 2, MPI_DOUBLE, _neighbours[1], 3,
			MPI_COMM_WORLD, &high);

	MPI_Isend(_leftOccBuffer, _xHaloOccSize, MPI_INT, _neighbours[0], 2,
				MPI_COMM_WORLD, &lowOcc);
	MPI_Isend(_rightOccBuffer, _xHaloOccSize, MPI_INT, _neighbours[1], 4,
				MPI_COMM_WORLD, &highOcc);

	MPI_Recv(_leftBufferRec, _xHaloSize * 2,MPI_DOUBLE, _neighbours[0],3,MPI_COMM_WORLD, &lowRecv);
	MPI_Recv(_rightBufferRec, _xHaloSize * 2,MPI_DOUBLE, _neighbours[1],1,MPI_COMM_WORLD, &highRecv);

	MPI_Recv(_leftOccBufferRec, _xHaloOccSize,MPI_INT, _neighbours[0],4,MPI_COMM_WORLD, &lowOccRecv);
	MPI_Recv(_rightOccBufferRec, _xHaloOccSize,MPI_INT, _neighbours[1],2,MPI_COMM_WORLD, &highOccRecv);

#endif
}

void UniformPseudoParticleContainer::communicateHalosY(){
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	MPI_Request low, high;
	MPI_Request lowOcc, highOcc;

	MPI_Status lowRecv,highRecv;
	MPI_Status lowOccRecv,highOccRecv;

	MPI_Isend(_bottomBuffer, _yHaloSize * 2, MPI_DOUBLE, _neighbours[2], 1,
			MPI_COMM_WORLD, &low);
	MPI_Isend(_topBuffer, _yHaloSize * 2, MPI_DOUBLE, _neighbours[3], 3,
			MPI_COMM_WORLD, &high);

	MPI_Isend(_bottomOccBuffer, _yHaloOccSize, MPI_INT, _neighbours[2], 2,
			MPI_COMM_WORLD, &lowOcc);
	MPI_Isend(_topOccBuffer, _yHaloOccSize, MPI_INT, _neighbours[3], 4,
			MPI_COMM_WORLD, &highOcc);

	MPI_Recv(_bottomBufferRec, _yHaloSize * 2, MPI_DOUBLE, _neighbours[2], 3, MPI_COMM_WORLD, &lowRecv);
	MPI_Recv(_topBufferRec, _yHaloSize * 2, MPI_DOUBLE, _neighbours[3], 1, MPI_COMM_WORLD, &highRecv);

	MPI_Recv(_bottomOccBufferRec, _yHaloOccSize, MPI_INT, _neighbours[2], 4, MPI_COMM_WORLD, &lowOccRecv);
	MPI_Recv(_topOccBufferRec, _yHaloOccSize, MPI_INT, _neighbours[3], 2, MPI_COMM_WORLD, &highOccRecv);

#endif
}

void UniformPseudoParticleContainer::communicateHalosZ(){
#if defined(ENABLE_MPI) && defined(NEW_FMM)
	MPI_Request low, high;
	MPI_Request lowOcc, highOcc;

	MPI_Status lowRecv,highRecv;
	MPI_Status lowOccRecv,highOccRecv;

	MPI_Isend(_backBuffer, _zHaloSize * 2, MPI_DOUBLE, _neighbours[4], 1,
			MPI_COMM_WORLD, &low);
	MPI_Isend(_frontBuffer, _zHaloSize * 2, MPI_DOUBLE, _neighbours[5], 3,
			MPI_COMM_WORLD, &high);

	MPI_Isend(_backOccBuffer, _zHaloOccSize, MPI_INT, _neighbours[4], 2,
			MPI_COMM_WORLD, &lowOcc);
	MPI_Isend(_frontOccBuffer, _zHaloOccSize, MPI_INT, _neighbours[5], 4,
			MPI_COMM_WORLD, &highOcc);

	MPI_Recv(_backBufferRec, _zHaloSize * 2, MPI_DOUBLE, _neighbours[4], 3, MPI_COMM_WORLD, &lowRecv);
	MPI_Recv(_frontBufferRec, _zHaloSize * 2, MPI_DOUBLE, _neighbours[5], 1, MPI_COMM_WORLD, &highRecv);

	MPI_Recv(_backOccBufferRec, _zHaloOccSize, MPI_INT, _neighbours[4], 4, MPI_COMM_WORLD, &lowOccRecv);
	MPI_Recv(_frontOccBufferRec, _zHaloOccSize, MPI_INT, _neighbours[5], 2, MPI_COMM_WORLD, &highOccRecv);

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
#if WIGNER==0
		GatherWellSepLo(cellWid, curCellsEdge, curLevel);
#else
		GatherWellSepLo_Wigner(cellWid, curCellsEdge, curLevel);
#endif
		AllReduceLocalMoments(curCellsEdge, curLevel);

		if(curLevel<_maxLevel) {
#if WIGNER==0
			PropagateCellLo(cellWid, curCellsEdge, curLevel);
#else
			PropagateCellLo_Wigner(cellWid, curCellsEdge, curLevel);
#endif
		}
	}
}

void UniformPseudoParticleContainer::processTreeInitM2LWigner() {


	int curCellsEdge = 1;
	double cellWid[3];

	for(int i=0; i <3; i++) cellWid[i] = _domain->getGlobalLength(i);

	for(int curLevel=1; curLevel<=_maxLevel; curLevel++){
		curCellsEdge *=2;
		for(int i=0; i <3; i++){
			cellWid[i] /=2;
		}
		GatherWellSepLoInitM2LWigner(cellWid, curCellsEdge, curLevel);

	}
}

void UniformPseudoParticleContainer::GatherWellSepLoInitM2LWigner(double *cellWid, int mpCells, int& curLevel) {
	int m1v[3];
	int m2v[3];

	int m1, m2, m2x, m2y, m2z;
	// int m1x, m1y, m1z;
	int m22x, m22y, m22z; // for periodic image
	int _size, _rank, loop_min, loop_max;
	int _row_length;

	_row_length=mpCells*mpCells*mpCells;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	_rank= domainDecomp.getRank();
	_size= domainDecomp.getNumProcs();

	loop_min = (int) ((long) (_rank + 0) * (long) (_row_length) / (long) _size);
	loop_max = (int) ((long) (_rank + 1) * (long) (_row_length) / (long) _size);

	Vector3<double> periodicShift;

	for (m1 = loop_min; m1 < loop_max; m1++) {

		m1v[0] = m1 % mpCells;
		m1v[1] = (m1 / mpCells) % mpCells;
		m1v[2] = (m1 / (mpCells * mpCells)) % mpCells;
//		if (_mpCell[curLevel][m1].occ == 0)
//			continue;

		for (m2z = LoLim(2); m2z <= HiLim(2); m2z++) {
			// to get periodic image
			m22z = (mpCells + m2z) % mpCells;
			periodicShift[2] = 0.0;
			if (m2z < 0) 		periodicShift[2] = -mpCells * cellWid[2];
			if (m2z >= mpCells) periodicShift[2] = mpCells * cellWid[2];

			m2v[2] = m2z;
			for (m2y = LoLim(1); m2y <= HiLim(1); m2y++) {
				// to get periodic image
				m22y = (mpCells + m2y) % mpCells;

				periodicShift[1] = 0.0;
				if (m2y < 0)		periodicShift[1] = -mpCells * cellWid[1];
				if (m2y >= mpCells)	periodicShift[1] = mpCells * cellWid[1];

				m2v[1] = m2y;
				for (m2x = LoLim(0); m2x <= HiLim(0); m2x++) {
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

					const SHMultipoleParticle& sh_multipole = _mpCell[curLevel][m2].multipole;
					const SHLocalParticle& sh_local = _mpCell[curLevel][m1].local;

					// compute periodically-shifted center
					Vector3<double> shifted_center = sh_multipole.getCenter() + periodicShift;

					// distance-vector FROM local TO multipole
					Vector3<double> r_target_to_source = shifted_center - sh_local.getCenter();

					/* compute index vector and add param to map */
					Vector3<int> idxVec(rint(r_target_to_source[0]/(cellWid[0])),
							rint(r_target_to_source[1]/(cellWid[1])),
							rint(r_target_to_source[2]/(cellWid[2])));

					// continue if element is already in container
					if (_M2L_Wigner.find(idxVec) != _M2L_Wigner.end() ) continue;

					double projXY_len = sqrt(
							r_target_to_source[0] * r_target_to_source[0]
									+ r_target_to_source[1] * r_target_to_source[1]);

					double phi = atan2(r_target_to_source[1], r_target_to_source[0]);
					double theta = atan2(projXY_len, r_target_to_source[2]);

					double* CosSin = new double[(_maxOrd+1)*2];
					for (int m = 0; m <= _maxOrd; ++m) {
						CosSin[2*m] 		= cos(m*phi);
						CosSin[2*m + 1] 	= sin(m*phi);
					}

					WignerMatrix W_pos(_maxOrd, true);
					WignerMatrix W_neg(_maxOrd, true);

					W_pos.evaluate(theta);
					W_neg.evaluate(-theta);

					// pre-multiply prefactors
					for (int l = 0; l <= _maxOrd; ++l) {
						for (int m = 0; m <= l; ++m) {
							for (int k = -l; k<=l; ++k) {
								W_pos.acc(l,m,k) *= bhfmm::RotationParameterLookUp::tab->acc_c(l,m,k);
								W_neg.acc(l,m,k) *= bhfmm::RotationParameterLookUp::tab->acc_c(l,k,m);
							}
						}
					}

					RotationParams& param = _M2L_Wigner[idxVec];// {W_pos, W_neg, CosSin};
					param.W[0] = W_pos;
					param.W[1] = W_neg;
					param.SinCos = CosSin;

				} // m2x closed
			} // m2y closed
		} // m2z closed
	} //m1 closed


}

void UniformPseudoParticleContainer::printTimers() {
	std::cout << "\t\t" << _timerAllreduce.get_etime()       		<< "\t\t" <<"s in Allreduce" << std::endl;
	std::cout << "\t\t" << _timerAllreduce_me.get_etime()			<< "\t\t" <<"s in Allreduce_me"<<std::endl;
	std::cout << "\t\t" << _timerCombineMpCell.get_etime()     		<< "\t\t" <<"s in CombineMpCell" << std::endl;
	std::cout << "\t\t" << _timerGatherWellSepLo.get_etime() 		<< "\t\t" <<"s in GatherWellSepLo" << std::endl;
	std::cout << "\t\t" << _timerPropagateCellLo.get_etime() 		<< "\t\t" <<"s in PropagateCellLo" << std::endl;
}


} /* namespace bhfmm */


