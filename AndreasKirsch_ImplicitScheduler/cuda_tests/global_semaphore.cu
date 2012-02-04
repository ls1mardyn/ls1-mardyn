#include <cuda_runtime.h>
#include <stdio.h>

class Mutex {
public:
	volatile uint lockValue;

public:
	__device__ void init() {
		lockValue = 0u;
	}

	// locks but supports being interrupted when (signal & bitMask) becomes 0
	// returns true if unlocked must be called
	// returns false if unlocked mustn't be called
	// assumption:
	// the signal can only be reset by a warp/thread that is inside the same semaphore,
	// only the caller can set it
	__device__ bool lock(volatile uint *signal, uint bitMask) {
		while( atomicExch( (uint*) &lockValue, 1u ) != 0u ) {
			if( (*signal & bitMask) == 0 ) {
				return false;
			}
		}

		return true;
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

volatile __shared__ uint warpSignals;

static __global__ void kernel() {
	if( threadIdx.x + threadIdx.y == 0 ) {
		printf( "warp 1 lag\n" );
	}

	if( threadIdx.x + threadIdx.y == 0 ) {
		warpSignals = 0u;
	}
	__syncthreads();
	// the __syncthreads and the warp == 1 init is necessary otherwise warpSignals might be reset later on
	// or will not have been initialized properly when it is accessed in other threads

	const uint bitMask = 1 << threadIdx.y;

	if( threadIdx.x == 0 ) {
		atomicOr( (uint*) &warpSignals, bitMask );

		//WARP_PRINTF( "signal: %i\n", warpSignals );
		//WARP_PRINTF( "signal: %i\n", warpSignals );

		if( globalMutex.lock( &warpSignals, bitMask ) ) {
			//WARP_PRINTF( "got lock\n" );

			uint signalSnapshot = warpSignals;

			globalCounter += __popc(signalSnapshot);

			//WARP_PRINTF( "signals\n" );
			atomicXor( (uint*) &warpSignals, signalSnapshot );

			//WARP_PRINTF( "release lock\n" );
			globalMutex.unlock();
		}
		else {
			//WARP_PRINTF( "got signal\n" );
		}
	}
}

void testGlobalSemaphore() {
	init<<<1,1>>>();
	dim3 blockDim = dim3( 32, 16, 1 );
	kernel<<<16, blockDim>>>();
	printResults<<<1,1>>>();
}
