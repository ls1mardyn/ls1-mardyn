// Andreas Kirsch <kirschan@tum.de>

/*
 * LinkedCellsOpenCL.h
 *
 *  Created on: Jul 4, 2010
 *      Author: orend
 */

#ifndef LINKEDCELLSOPENCL_H_
#define LINKEDCELLSOPENCL_H_
#define __NO_STD_VECTOR
#define __CL_ENABLE_EXCEPTIONS
//#define DISTANCES
//#define DEBUG_DISTANCES
#include "LinkedCells.h"
#include "Cell.h"
#include "ParticleContainer.h"
#include <malloc.h>

class LinkedCellsOpenCL: public LinkedCells {
public:
	LinkedCellsOpenCL(double bBoxMin[3], double bBoxMax[3], double cutoffRadius, double LJCutoffRadius,
		     double tersoffCutoffRadius, double cellsInCutoffRadius,
		     ParticlePairsHandler* partPairsHandler);
	virtual ~LinkedCellsOpenCL();
	void traversePairs();

private:
	// no CUDA stream support (yet?)

	// host memory
	int* pairs;
	float* m_positions;
	float* distances;
	float* forces;

	int numberOfPairs;
	int numberOfSelfs;
	int resultSize;
	int seflResultSize;
	//int maxCellSize;
	const static int max_threads=512;
	double cutoffRadiusSquare;
	double LJCutoffRadiusSquare;
	double tersoffCutoffRadiusSquare;
	double totalTime;
	long numberOfParticles;
	void initCUDA();
	void traversePairsInit();
	void traversePairsFinish();
	int getMaxCellSizeAndPositions();
	int countParticles();
	void countPairs();
	void createPairs();
	void invertPairs();
};

#endif /* LINKEDCELLSOPENCL_H_ */
