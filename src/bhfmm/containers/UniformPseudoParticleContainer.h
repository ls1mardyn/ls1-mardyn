/*
 * UniformPseudoParticleContainer.h
 *
 *  Created on: Feb 5, 2015
 *      Author: tchipevn
 */

#ifndef UNIFORMPSEUDOPARTICLECONTAINER_H_
#define UNIFORMPSEUDOPARTICLECONTAINER_H_

#include "PseudoParticleContainer.h"
#include "LeafNodesContainer.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Timer.h"
#include "bhfmm/utils/WignerMatrix.h"
#include "bhfmm/utils/RotationParameter.h"
#include <vector>
#include <map>
#include "bhfmm/HaloBufferNoOverlap.h"
#include "bhfmm/HaloBufferOverlap.h"

class Domain;
class DomainDecompBase;

namespace bhfmm {

class UniformPseudoParticleContainer: public PseudoParticleContainer {
public:
	UniformPseudoParticleContainer(double domainLength[3], double bBoxMin[3], double bBoxMax[3],
			double LJCellLength[3], unsigned LJSubdivisionFactor, int orderOfExpansions,bool periodic = true);
	~UniformPseudoParticleContainer();

	void clear();
	void build(ParticleContainer* pc);
	void upwardPass(P2MCellProcessor * cp);
	void horizontalPass(VectorizedChargeP2PCellProcessor * cp);
	void downwardPass(L2PCellProcessor *cp);

	// P2M
	void processMultipole(ParticleCell& cell);

	// L2P
	void processFarField(ParticleCell& cell);

	// M2M, M2L, L2L
	void processTree();

	void printTimers();



private:
	LeafNodesContainer* _leafContainer;

	int _wellSep;
	int _maxLevel;
	int _globalLevel;
	//In the parallel version the octree is divided into two trees:
	//- A local subtree starting at _globalLevel + 1 which contains only the
	//ancestor nodes from the node that was assigned to the MPI process on
	//the _globalLevel
	//- A global tree which contains all the octree elements of the FMM tree
	//from level 0 to _globalLevel
	std::vector<std::vector<MpCell> > _mpCellGlobalTop;
	std::vector<std::vector<MpCell> > _mpCellLocal;
	MPI_Request _allReduceRequest;
	double _cellLength[3];
	int _globalNumCellsPerDim;
	Domain* _domain;
	int _globalNumCells;
	int* _occVector; // array for MPI allgather

	int _coeffVectorLength;
	double* _coeffVector; // array for MPI allgather
	double* _coeffVector_me;

	HaloBufferNoOverlap<double> * _multipoleRecBuffer, *_multipoleBuffer;
	HaloBufferOverlap<double> * _multipoleRecBufferOverlap, *_multipoleBufferOverlap;

	bool _periodicBC;

	int _numProcessorsPerDim;
	Vector3<int> _processorPositionGlobalLevel;
	Vector3<double> _bBoxMin;
	std::vector<int> _neighbours;

	int _globalLevelNumCells;
	// M2M
	void CombineMpCell(double *cellWid, int& mpCells, int curLevel);

	// M2M
	void CombineMpCell_MPI(double *cellWid, int& mpCells, int curLevel, Vector3<int> offset);

	// M2L
	void GatherWellSepLo(double *cellWid, int mpCells, int curLevel);

	// M2L
	void GatherWellSepLo_MPI(double *cellWid, int mpCells, int curLevel, int doHalos);

	// L2L
	void PropagateCellLo(double *cellWid, int mpCells, int curLevel);

	// L2L
	void PropagateCellLo_MPI(double *cellWid, int mpCells, int curLevel, Vector3<int> offset);

	// for parallelization
	void AllReduceMultipoleMoments();
	void AllReduceLocalMoments(int mpCells, int _curLevel);
	void AllReduceMultipoleMomentsLevelToTop(int mpCells, int _curLevel);
	void AllReduceMultipoleMomentsSetValues(int mpCells, int _curLevel);

	void getXHaloValues(int localMpCellsBottom,int bottomLevel);
	void getYHaloValues(int localMpCellsBottom,int bottomLevel);
	void getZHaloValues(int localMpCellsBottom,int bottomLevel);
	void setXHaloValues(int localMpCellsBottom,int bottomLevel);
	void setYHaloValues(int localMpCellsBottom,int bottomLevel);
	void setZHaloValues(int localMpCellsBottom,int bottomLevel);
	void getHaloValues(int localMpCellsBottom,int bottomLevel, double *buffer,
			int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh);
	void setHaloValues(int localMpCellsBottom,int bottomLevel, double *bufferRec,
			int xLow, int xHigh, int yLow, int yHigh, int zLow, int zHigh);
	//for parallelization
	void communicateHalosNoOverlap();
	void communicateHalosOverlapStart();
	//has to be called after receive finished
	void communicateHalosOverlapSetHalos();
	void communicateHalos();
	void communicateHalosX();
	void communicateHalosY();
	void communicateHalosZ();
	void communicateHalosAlongAxis(double * lowerNeighbourBuffer, double * higherNeighbourBuffer,
			double * lowerNeighbourBufferRec, double * higherNeighbourBufferRec,
			int lowerNeighbour, int higherNeighbour, int haloSize
			);


	Timer _timerProcessCells;
	Timer _timerAllreduce;
	Timer _timerCombineMpCell;
	Timer _timerGatherWellSepLo;
	Timer _timerPropagateCellLo;
	Timer _timerProcessFarField;

	Timer _timerGatherEvalM;
	Timer _timerGatherEvalLM;
	Timer _timerAllreduce_me;

	MPI_Comm _comm;
	int _overlapComm;

};

} /* namespace bhfmm */

#endif /* UNIFORMPSEUDOPARTICLECONTAINER_H_ */
