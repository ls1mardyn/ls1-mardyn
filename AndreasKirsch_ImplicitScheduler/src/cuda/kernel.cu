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

// use these printf wrappers for simplified debug output
#if 0
#	define THREAD_PRINTF(format, ...) printf( "({%i} %i W%i T%i) " format, blockIdx.x + blockIdx.y * gridDim.x, warpIdx, warpThreadIdx, getSM(), ##__VA_ARGS__ )
#	define WARP_PRINTF(format, ...) do { if( warpThreadIdx == 0 ) printf( "({%i} %i W%i) " format, blockIdx.x + blockIdx.y * gridDim.x, warpIdx, getSM(), ##__VA_ARGS__ ); } while( false )
#	define BLOCK_PRINTF(format, ...) do { if( warpIdx == 0 ) printf( "({%i} %i) " format, blockIdx.x + blockIdx.y * gridDim.x, getSM(), ##__VA_ARGS__ ); } while( false )
#	define GRID_PRINTF(format, ...) do { if( blockIdx.x == 0 && blockIdx.y == 0 && warpThreadIdx == 0 && warpIdx == 0 ) printf( "{%i} " format, getSM(), ##__VA_ARGS__ ); } while( false )
#else
#	define THREAD_PRINTF(format, ...)
#	define WARP_PRINTF(format, ...)
#	define BLOCK_PRINTF(format, ...)
#	define GRID_PRINTF(format, ...)
#endif


#include "domainTraverser.cum"
#include "globalStats.cum"

#include "referenceCellProcessor.cum"
#include "threadBlockCellProcessor.cum"
#include "warpBlockCellProcessor.cum"

#ifdef REFERENCE_IMPLEMENTATION
#warning using reference cell processor
#endif

#ifdef CUDA_DOUBLE_MODE
#	warning using double precision
#else
#	warning using float precision
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

	const Quaternion quaternion = moleculeQuaternions[ moleculeIndex ];

#ifndef TEST_QUATERNION_MATRIX_CONVERSION
	moleculeRotations[ moleculeIndex ] = quaternion.toRotMatrix3x3();
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

#ifndef REFERENCE_IMPLEMENTATION
namespace CellProcessor = ThreadBlockCellProcessor;
#else
namespace CellProcessor = ReferenceCellProcessor;
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
	CellProcessor::processCellPair( threadIndex, cellA, cellB );

	ThreadBlockCellStats::reduceAndStore( threadIndex, cellIndex, neighborIndex );
}

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
	CellProcessor::processCell( threadIndex, cell );

	ThreadBlockCellStats::reduceAndStore( threadIndex, cellIndex, cellIndex );
}
#else
namespace CellProcessor = WarpBlockCellProcessor;

__device__ CellProcessor::CellScheduler *cellScheduler;
__device__ CellProcessor::CellPairScheduler *cellPairScheduler;

__global__ void createSchedulers() {
	cellScheduler = new CellProcessor::CellScheduler();
	cellPairScheduler = new CellProcessor::CellPairScheduler();
}

__global__ void destroySchedulers() {
	delete cellScheduler;
	delete cellPairScheduler;
}

__global__ void processCellPair() {
	const int threadIndex = getThreadIndex();

	__shared__ CellProcessor::ThreadBlockInfo threadBlockInfo;
	if( threadIndex == 0 ) {
		threadBlockInfo.init();
	}
	__syncthreads();

	do {
		cellPairScheduler->scheduleWarpBlocks( threadBlockInfo );

		while( !threadBlockInfo.warpJobQueue[warpIdx].isEmpty() ) {
			ThreadBlockCellStats::initThreadLocal( threadIndex );

			CellProcessor::WarpBlockPairInfo warpBlockPairInfo = threadBlockInfo.warpJobQueue[warpIdx].pop();
			CellProcessor::processCellPair( warpBlockPairInfo );

			ThreadBlockCellStats::reduceAndStoreWarp( threadIndex, warpBlockPairInfo.warpBlockA.cellIndex );
		}
	} while( threadBlockInfo.hasMoreJobs );

	WARP_PRINTF( "terminating\n" );
}

__global__ void processCell() {
	// TODO: remove?
	const int threadIndex = getThreadIndex();

	__shared__ CellProcessor::ThreadBlockInfo threadBlockInfo;
	if( threadIndex == 0 ) {
		threadBlockInfo.init();
	}
	__syncthreads();

	do {
		cellScheduler->scheduleWarpBlocks( threadBlockInfo );

		while( !threadBlockInfo.warpJobQueue[warpIdx].isEmpty() ) {
			ThreadBlockCellStats::initThreadLocal( threadIndex );

			CellProcessor::WarpBlockPairInfo warpBlockPairInfo = threadBlockInfo.warpJobQueue[warpIdx].pop();
			CellProcessor::processCell( warpBlockPairInfo );

			ThreadBlockCellStats::reduceAndStoreWarp( threadIndex, warpBlockPairInfo.warpBlockA.cellIndex );
		}
	} while( threadBlockInfo.hasMoreJobs );

	WARP_PRINTF( "terminating\n" );
}

#endif

}
