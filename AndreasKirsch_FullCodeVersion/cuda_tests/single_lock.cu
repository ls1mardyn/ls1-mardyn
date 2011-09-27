#include <cuda_runtime.h>
#include <stdio.h>

class Mutex {
public:
	volatile uint lockValue;

public:
	__device__ void init() {
		lockValue = 0u;
	}

	__device__ void lock() {
		while( atomicExch( (uint*) &lockValue, 1u ) != 0u )
			;
	}

	__device__ void unlock() {
		lockValue = 0u;
	}
};

static __device__ Mutex globalMutex;
static volatile __device__ uint globalCounter;

#define WARP_PRINTF( fmt, ... ) printf( "(%i %i) " fmt, blockIdx.x, threadIdx.y, ##__VA_ARGS__ )

static __global__ void init() {
	globalMutex.lockValue = 0;
	globalCounter = 0;
}

static __global__ void printResults() {
	printf( "%i\n", globalCounter );
}

static __global__ void kernel() {
	if( threadIdx.x == 0 ) {
		globalMutex.lock();

		//WARP_PRINTF( "got lock\n" );
		globalCounter++;

		//WARP_PRINTF( "incremented globalCounter\n" );
		__threadfence();

		//WARP_PRINTF( "releasing lock\n" );
		globalMutex.unlock();
	}
}

void testSingleLock() {
	init<<<1,1>>>();
	dim3 blockDim = dim3( 32, 16, 1 );
	kernel<<<16, blockDim>>>();
	printResults<<<1,1>>>();
}

