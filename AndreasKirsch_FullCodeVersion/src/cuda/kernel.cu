#include <host_defines.h>
#include <stdio.h>

#include "cutil_math.h"

#include "config.h"

__device__ int getSM() {
	uint smID;

	asm volatile( "mov.u32 %0, %smid;" : "=r"(smID) );

	return smID;
}

#define warpThreadIdx threadIdx.x
#define warpIdx threadIdx.y

__device__ __forceinline__ uint getThreadIndex() {
	return warpThreadIdx + WARP_SIZE * warpIdx;
}

#if 0
#	define WARP_PRINTF(format, ...) if( warpThreadIdx == 0 ) printf( "(%i W%i {%i}) " format, blockIdx.x + blockIdx.y * gridDim.x, warpIdx, getSM(), ##__VA_ARGS__ )
#	define BLOCK_PRINTF(format, ...) if( warpIdx == 0 ) printf( "(%i {%i}) " format, blockIdx.x + blockIdx.y * gridDim.x, getSM(), ##__VA_ARGS__ )
#	define GRID_PRINTF(format, ...) if( blockIdx.x == 0 && blockIdx.y == 0 && warpThreadIdx == 0 && warpIdx == 0 ) printf( "{%i} " format, getSM(), ##__VA_ARGS__ )
#else
#	define WARP_PRINTF(format, ...)
#	define BLOCK_PRINTF(format, ...)
#	define GRID_PRINTF(format, ...)
#endif

#include "moleculeStorage.cum"

#include "componentDescriptor.cum"

#include "moleculePairHandler.cum"

#include "domainTraverser.cum"

#include "referenceCellProcessor.cum"

#include "threadBlockCellProcessor.cum"

#include "globalStats.cum"

#include "cellInfo.cum"

#include "molecule.cum"

#include "warpBlockCellProcessor.cum"

#ifndef REFERENCE_IMPLEMENTATION
#warning using fast cell processor
#else
#warning using reference cell processor
#endif

#ifdef CUDA_DOUBLE_MODE
#	warning using double precision
#else
#	warning using float precision
#endif

#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
#	warning sorting cells by component type
#else
#	warning cells are *not* sorted by component type
#endif

#ifdef CUDA_HW_CACHE_ONLY
#	warning no shared local storage cache
#else
#	warning shared local storage active
#endif

#ifdef NO_CONSTANT_MEMORY
#	warning no constant memory
#else
#	warning constant memory used
#endif

extern "C" {
/* TODO: possible refactoring
 * create a prepare method in MoleculeStorage
 */
// TODO: interesting to benchmark idea: unrolled loop in this kernel vs the way it is now---overhead?
__global__ void convertQuaternionsToRotations( int numMolecules ) {
	int moleculeIndex = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	if( moleculeIndex >= numMolecules ) {
		return;
	}

#ifndef CUDA_UNPACKED_STORAGE
	const Quaternion quaternion = moleculeQuaternions[ moleculeIndex ];
#else
	const Quaternion quaternion = packQuaternion( moleculeQuaternions, moleculeIndex );
#endif

#ifndef TEST_QUATERNION_MATRIX_CONVERSION
#	ifndef CUDA_UNPACKED_STORAGE
	moleculeRotations[ moleculeIndex ] = quaternion.toRotMatrix3x3();
#	else
	unpackMatrix3x3( moleculeRotations, moleculeIndex, quaternion.toRotMatrix3x3() );
#	endif
#else
#	warning CUDA: testing quaternion matrix conversion
	const Matrix3x3 convertedQuaternion = quaternion.toRotMatrix3x3();
	const Matrix3x3 &correctRotation = moleculeRotations[ moleculeIndex ];

	const floatType error = length( convertedQuaternion.rows[0] - correctRotation.rows[0] ) +
			length( convertedQuaternion.rows[1] - correctRotation.rows[1] ) +
			length( convertedQuaternion.rows[2] - correctRotation.rows[2] );

	if( error > 1e-9 ) {
		printf( "bad quaternion conversion (molecule %i)\n", moleculeIndex );
	}
#endif
}

#ifndef CUDA_WARP_BLOCK_CELL_PROCESSOR

__device__ MoleculeStorage moleculeStorage;

#ifndef REFERENCE_IMPLEMENTATION
#	ifndef CUDA_HW_CACHE_ONLY
__shared__ SharedMoleculeLocalStorage< moleculeStorage > moleculeLocalStorage;
#	else
__device__ WriteThroughMoleculeLocalStorage< moleculeStorage > moleculeLocalStorage;
#	endif
__device__ ThreadBlockCellProcessor<
	typeof(moleculeStorage), moleculeStorage,
	typeof(moleculeLocalStorage), moleculeLocalStorage>
		cellProcessor;

#else
__device__ ReferenceCellProcessor<
	typeof(moleculeStorage), moleculeStorage,
	typeof(moleculePairHandler), moleculePairHandler>
		cellProcessor;
#endif

__global__ void processCellPair() {
	const int threadIndex = getThreadIndex();

	const int jobIndex = blockIdx.y * gridDim.x + blockIdx.x;
	if( jobIndex >= DomainTraverser::numJobs ) {
		return;
	}

	int cellIndex = DomainTraverser::getCellIndexFromJobIndex( jobIndex );
	int neighborIndex = DomainTraverser::getNeighborCellIndex( cellIndex );

	CellInfoEx cellA = cellInfoFromCellIndex( cellIndex );
	CellInfoEx cellB = cellInfoFromCellIndex( neighborIndex );

	if( cellA.length == 0 || cellB.length == 0 ) {
		return;
	}

	ThreadBlockCellStats::initThreadLocal( threadIndex );
	cellProcessor.processCellPair( threadIndex, cellA, cellB );

	ThreadBlockCellStats::reduceAndStore( threadIndex, cellIndex, neighborIndex );
}

//__launch_bounds__(BLOCK_SIZE, 2)
__global__ void processCell() {
	const int threadIndex = getThreadIndex();

	int jobIndex = blockIdx.y * gridDim.x + blockIdx.x;
	if( jobIndex >= DomainTraverser::numJobs ) {
		return;
	}

	int cellIndex = DomainTraverser::getInnerCellIndexFromJobIndex(jobIndex);
	CellInfoEx cell = cellInfoFromCellIndex( cellIndex );
	if( cell.length == 0 ) {
		return;
	}

	ThreadBlockCellStats::initThreadLocal( threadIndex );
	cellProcessor.processCell( threadIndex, cell );

	ThreadBlockCellStats::reduceAndStore( threadIndex, cellIndex, cellIndex );
}
#else

__device__ WBCP::CellScheduler *cellScheduler;
__device__ WBCP::CellPairScheduler *cellPairScheduler;

__global__ void createSchedulers() {
	cellScheduler = new WBDP::CellScheduler();
	cellPairScheduler = new WBDP::CellPairScheduler();
}

__global__ void destroySchedulers() {
	delete cellScheduler;
	delete cellPairScheduler;
}

__device__ MoleculeStorage moleculeStorage;
__shared__ ResultLocalStorage< moleculeStorage > resultLocalStorage;
__device__ WBCP::CellProcessor< moleculeStorage, typeof(resultLocalStorage), resultLocalStorage > cellProcessor;

__global__ void processCellPair() {
	const int threadIndex = getThreadIndex();

	__shared__ WBCP::ThreadBlockInfo threadBlockInfo;
	if( threadIndex == 0 ) {
		threadBlockInfo.init();
	}
	__syncthreads();

	do {
		cellPairScheduler->scheduleWarpBlocks( threadBlockInfo );

		while( !threadBlockInfo.warpJobQueue[warpIdx].isEmpty() ) {
			ThreadBlockCellStats::initThreadLocal( threadIndex );

			WBCP::WarpBlockPairInfo warpBlockPairInfo = threadBlockInfo.warpJobQueue[warpIdx].pop();
#ifdef CUDA_HW_CACHE_ONLY
			cellProcessor.processCellPair( warpBlockPairInfo );

			ThreadBlockCellStats::reduceAndStoreWarp( threadIndex, warpBlockPairInfo.warpBlockA.cellIndex );
#else
			cellProcessor.processCellPairWithCache( warpBlockPairInfo );

			ThreadBlockCellStats::reduceAndStoreWarpForPair( threadIndex, warpBlockPairInfo.warpBlockA.cellIndex );
#endif
		}
	} while( threadBlockInfo.hasMoreJobs );

	WARP_PRINTF( "terminating..\n" );
}

__global__ void processCell() {
	// TODO: remove?
	const int threadIndex = getThreadIndex();

	__shared__ WBCP::ThreadBlockInfo threadBlockInfo;
	if( threadIndex == 0 ) {
		threadBlockInfo.init();
	}
	__syncthreads();

	do {
		cellScheduler->scheduleWarpBlocks( threadBlockInfo );

		while( !threadBlockInfo.warpJobQueue[warpIdx].isEmpty() ) {
			ThreadBlockCellStats::initThreadLocal( threadIndex );

			WBCP::WarpBlockPairInfo warpBlockPairInfo = threadBlockInfo.warpJobQueue[warpIdx].pop();
			cellProcessor.processCell( warpBlockPairInfo );

			ThreadBlockCellStats::reduceAndStoreWarp( threadIndex, warpBlockPairInfo.warpBlockA.cellIndex );
		}
	} while( threadBlockInfo.hasMoreJobs );

	WARP_PRINTF( "terminating\n" );
}

#endif

}
