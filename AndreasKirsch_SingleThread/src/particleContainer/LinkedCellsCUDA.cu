// Andreas Kirsch 2010

#include "LinkedCellsCUDA.h"
#include "molecules/Molecule.h"
#include "cutil_math.h"

#define OUT

#define CUDA_TIMING

#ifdef CUDA_TIMING
class CUDATimer {
private:
	cudaEvent_t _startEvent, _endEvent;

public:
	CUDATimer() {
		CUDA_THROW_ON_ERROR( cudaEventCreate( &_startEvent ) );
		CUDA_THROW_ON_ERROR( cudaEventCreate( &_endEvent ) );
	}

	~CUDATimer() {
		CUDA_THROW_ON_ERROR( cudaEventDestroy( _startEvent ) );
		CUDA_THROW_ON_ERROR( cudaEventDestroy( _endEvent ) );
	}

	void begin() {
		CUDA_THROW_ON_ERROR( cudaEventRecord( _startEvent ) );
	}

	void end() {
		CUDA_THROW_ON_ERROR( cudaEventRecord( _endEvent ) );
	}

	float getElapsedTime() {
		CUDA_THROW_ON_ERROR( cudaEventSynchronize( _endEvent ) );

		float elapsedTime;
		CUDA_THROW_ON_ERROR( cudaEventElapsedTime( &elapsedTime, _startEvent, _endEvent ) );

		return elapsedTime;
	}

	void printElapsedTime( const char *format ) {
		printf( format, getElapsedTime() );
	}
};
#else
class CUDATimer {
public:
	void begin() {
	}

	void end() {
	}

	float getElapsedTime() {
		return 0.0f;
	}

	void printElapsedTime( const char *format ) {
	}
};
#endif

__device__ void calculateLennardJones( const float3 distance, const float distanceSquared, float epsilon, float sigmaSquared,
		OUT float3 &force, OUT float &potential) {
	float invdr2 = 1.f / distanceSquared;
	float lj6 = sigmaSquared * invdr2; lj6 = lj6 * lj6 * lj6;
	float lj12 = lj6 * lj6;
	float lj12m6 = lj12 - lj6;
	potential = 4.0f * epsilon * lj12m6;
	// result: force = fac * distance = fac * |distance| * normalized(distance)
	float fac = -24.0f * epsilon * (lj12 + lj12m6) * invdr2;
	force = fac * distance;
}

__device__ int getCellIndex( int startIndex, int2 dimension, int3 gridOffsets ) {
	const int idx = blockIdx.x;
	const int3 gridIndex = make_int3(
			idx % dimension.x,
			(idx / dimension.x) % dimension.y,
			idx / dimension.x / dimension.y
		);
	const int cellIndex = startIndex + dot( gridIndex, gridOffsets );

	return cellIndex;
}

__global__ void Kernel_calculateInnerLJForces( float3 *positions, OUT float3 *forces, int2 *cellInfos, OUT float2 *domainValues,
		float epsilon, float sigmaSquared, float cutOffRadiusSquared ) {
	const int cellIndex = blockIdx.x;

	const int cellStart = cellInfos[cellIndex].x;
	const int cellLength = cellInfos[cellIndex].y;

	float totalPotential = 0.0f;
	float totalVirial = 0.0f;

	// start with 1 because aIndex = 1 will result in an empty loop cycle
	for( int aIndex = 1 ; aIndex < cellLength ; aIndex++ ) {
		for( int bIndex = 0 ; bIndex < aIndex ; bIndex++ ) {
			const float3 &aPosition = positions[cellStart + aIndex];
			const float3 &bPosition = positions[cellStart + bIndex];

			const float3 distance = bPosition - aPosition;
			const float distanceSquared = dot( distance, distance );
			if( distanceSquared > cutOffRadiusSquared ) {
				continue;
			}

			float3 force;
			float potential;
			calculateLennardJones( distance, distanceSquared, epsilon, sigmaSquared, force, potential );

			totalPotential += potential;
			float virial = dot( force, distance );
			totalVirial += virial;

			forces[cellStart + aIndex] += force;
			forces[cellStart + bIndex] -= force;
		}
	}
	domainValues[cellIndex].x = totalPotential * 2;
	domainValues[cellIndex].y = totalVirial * 2;
}

__global__ void Kernel_calculatePairLJForces( float3 *positions, OUT float3 *forces, int2 *cellInfos, OUT float2 *domainValues,
		int startIndex, int2 dimension, int3 gridOffsets,
		int neighborOffset,
		float epsilon, float sigmaSquared, float cutOffRadiusSquared ) {
	const int cellIndex = getCellIndex( startIndex, dimension, gridOffsets );
	const int neighborIndex = cellIndex + neighborOffset;

	const int cellA_start = cellInfos[cellIndex].x;
	const int cellA_length = cellInfos[cellIndex].y;
	const int cellB_start = cellInfos[neighborIndex].x;
	const int cellB_length = cellInfos[neighborIndex].y;

	float totalPotential = 0.0f;
	float totalVirial = 0.0f;

	for( int aIndex = 0 ; aIndex < cellA_length ; aIndex++ ) {
		for( int bIndex = 0 ; bIndex < cellB_length ; bIndex++ ) {
			const float3 &aPosition = positions[cellA_start + aIndex];
			const float3 &bPosition = positions[cellB_start + bIndex];

			const float3 distance = bPosition - aPosition;
			const float distanceSquared = dot( distance, distance );
			if( distanceSquared > cutOffRadiusSquared ) {
				continue;
			}

			float3 force;
			float potential;
			calculateLennardJones( distance, distanceSquared, epsilon, sigmaSquared, force, potential );

			totalPotential += potential;
			float virial = dot( force, distance );
			totalVirial += virial;

			forces[cellA_start + aIndex] += force;
			forces[cellB_start + bIndex] -= force;
		}
	}

	domainValues[cellIndex].x += totalPotential;
	domainValues[cellIndex].y += totalVirial;
	domainValues[neighborIndex].x += totalPotential;
	domainValues[neighborIndex].y += totalVirial;
}

LinkedCellsCUDA_Internal::DomainValues LinkedCellsCUDA_Internal::calculateForces() {
	manageAllocations();

	initCellInfosAndCopyPositions();
	prepareDeviceMemory();

	calculateAllLJFoces();

	DomainValues domainValues;
	extractResultsFromDeviceMemory();
	reducePotentialAndVirial( domainValues.potential, domainValues.virial );

	printf( "Potential: %f Virial: %f\n", domainValues.potential, domainValues.virial );
	printf( "Average Potential: %f Average Virial: %f\n", domainValues.potential / _numParticles, domainValues.virial / _numParticles );

	determineForceError();

	updateMoleculeForces();

	return domainValues;
}

void LinkedCellsCUDA_Internal::manageAllocations()
{
	_numParticles = _linkedCells.getParticles().size();
	_numCells = _linkedCells.getCells().size();

	// TODO: use memalign like the old code?
	if( _numParticles > _maxParticles ) {
		_positions.resize( _numParticles );
		_forces.resize( _numParticles );

		_maxParticles = _numParticles;
	}

	if( _numCells > _maxCells ) {
		_cellInfos.resize( _numCells );
		_domainValues.resize( _numCells );

		_maxCells = _numCells;
	}
}

void LinkedCellsCUDA_Internal::freeAllocations()
{
	_positions.resize( 0 );
	_forces.resize( 0 );

	_cellInfos.resize( 0 );
	_domainValues.resize( 0 );
}

void LinkedCellsCUDA_Internal::initCellInfosAndCopyPositions()
{
	int currentIndex = 0;
	for( int i = 0 ; i < _numCells ; i++ ) {
		const Cell &cell = _linkedCells.getCells()[i];

		_cellInfos[i].x = currentIndex;
		_cellInfos[i].y = cell.getMoleculeCount();

		const std::list<Molecule*> &particles = cell.getParticlePointers();
		for( std::list<Molecule*>::const_iterator iterator = particles.begin() ; iterator != particles.end() ; iterator++ ) {
			Molecule &molecule = **iterator;
			_positions[currentIndex].x = molecule.r(0);
			_positions[currentIndex].y = molecule.r(1);
			_positions[currentIndex].z = molecule.r(2);
			currentIndex++;
		}
	}
}

void LinkedCellsCUDA_Internal::prepareDeviceMemory()
{
	// TODO: use page-locked/mapped memory
	printf( "%i\n", _numParticles );

	CUDATimer copyTimer;

	copyTimer.begin();

	_positions.copyToDevice();
	_cellInfos.copyToDevice();

	// init device forces to 0
	_forces.zeroDevice();
	// not needed: _domainValues.zeroDevice();

	copyTimer.end();
	copyTimer.printElapsedTime( "host to device copying: %f ms\n" );
}

void LinkedCellsCUDA_Internal::extractResultsFromDeviceMemory() {
	_forces.copyToHost();
	_domainValues.copyToHost();
}

void LinkedCellsCUDA_Internal::updateMoleculeForces() {
	int currentIndex = 0;
	for( int i = 0 ; i < _numCells ; i++ ) {
		const Cell &cell = _linkedCells.getCells()[i];

		const std::list<Molecule*> &particles = cell.getParticlePointers();
		for( std::list<Molecule*>::const_iterator iterator = particles.begin() ; iterator != particles.end() ; iterator++ ) {
			Molecule &molecule = **iterator;
			molecule.Fljcenterset( 0, (float*) &_forces[currentIndex] );
			currentIndex++;
		}
	}
}

void LinkedCellsCUDA_Internal::determineForceError() {
	double totalError = 0.0;
	double totalRelativeError = 0.0;
	float epsilon = 5.96e-06f;

	int currentIndex = 0;
	for( int i = 0 ; i < _numCells ; i++ ) {
		const Cell &cell = _linkedCells.getCells()[i];

		const std::list<Molecule*> &particles = cell.getParticlePointers();
		for( std::list<Molecule*>::const_iterator iterator = particles.begin() ; iterator != particles.end() ; iterator++ ) {
			Molecule &molecule = **iterator;
			float3 &cudaForce = _forces[currentIndex];
			const double *cpuForceD = molecule.ljcenter_F(0);
			float3 cpuForce = make_float3( cpuForceD[0], cpuForceD[1], cpuForceD[2] );
			float3 deltaForce = cudaForce - cpuForce;

			float error = length( deltaForce );
			totalError += error;

			if( error > epsilon ) {
				float relativeError = error / length( cpuForce );
				totalRelativeError += relativeError;
			}

			currentIndex++;
		}
	}

	printf( "Average Error: %f\nAverage Relative Error: %f\n", totalError / _numParticles, totalRelativeError / _numParticles );
}

void LinkedCellsCUDA_Internal::calculateAllLJFoces() {
	CUDATimer singleCellsTimer, cellPairsTimer;

	// TODO: wtf? this is from the old code
	const float epsilon = 1.0f;
	const float sigmaSquared = 1.0f;
	const float cutOffRadiusSquared = _cutOffRadius * _cutOffRadius;

	singleCellsTimer.begin();

	// inner forces first
	Kernel_calculateInnerLJForces<<<_numCells, 1>>>(
			_positions.devicePtr(), _forces.devicePtr(),_cellInfos.devicePtr(), _domainValues.devicePtr(),
			epsilon, sigmaSquared, cutOffRadiusSquared
		);

	singleCellsTimer.end();

	// pair forces
	const int *dimensions = _linkedCells.getCellDimensions();
	assert( dimensions[0] >= 2 && dimensions[1] >= 2 && dimensions[2] >=2 );

	const int3 zero3 = {0,0,0};
	const int3 xDirection = {1,0,0};
	const int3 yDirection = {0,1,0};
	const int3 zDirection = {0,0,1};
	// always make sure that each direction contains one component == 1
	const int3 directions[] = {
			{1,0,0},{0,1,0},{0,0,1},
			{1,1,0},{1,0,1},{0,1,1},
			{-1,1,0},{-1,0,1},{0,-1,1},
			{-1,1,1},{1,-1,1},{1,1,-1},
			{1,1,1}
	};

	cellPairsTimer.begin();

	for( int i = 0 ; i < sizeof( directions ) / sizeof( directions[0] ) ; i++ ) {
		const int3 &direction = directions[i];
		// we are going to iterate over odd and even slices (either xy-, xz- or yz-slices)

		// define: the main direction is the normal of the slice plane

		int neighborOffset = getDirectionOffset( direction );

		// contains the oriented direction as if the main direction was (0,0,1)
		int3 localDirection;
		// dimensions as if the main direction was (0,0,1)
		int3 localDimensions;
		int3 gridOffsets;

		// determine the direction of the plane (xy, xz or yz)
		if( direction.x == 1 ) {
			// yz plane (main direction: x)
			localDirection = make_int3( direction.y, direction.z, direction.x );
			localDimensions = make_int3( dimensions[1], dimensions[2], dimensions[0] );
			gridOffsets = make_int3(
					getDirectionOffset( yDirection ),
					getDirectionOffset( zDirection ),
					getDirectionOffset( xDirection )
				);
		}
		else if( direction.y == 1 ) {
			// xz plane (main direction: y)
			localDirection = make_int3( direction.x, direction.z, direction.y );
			localDimensions = make_int3( dimensions[0], dimensions[2], dimensions[1] );
			gridOffsets = make_int3(
					getDirectionOffset( xDirection ),
					getDirectionOffset( zDirection ),
					getDirectionOffset( yDirection )
				);
		}
		else if( direction.z == 1 ) {
			// xy plane (main direction: z)
			localDirection = direction;
			localDimensions = make_int3( dimensions[0], dimensions[1], dimensions[2] );
			gridOffsets = make_int3(
					getDirectionOffset( xDirection ),
					getDirectionOffset( yDirection ),
					getDirectionOffset( zDirection )
				);
		}
		else {
			assert( false );
		}

		// determine the startOffset as first cell near (0,0,0) so that start + neighborOffset won't be out of bounds
		int evenSlicesStartIndex = getCellOffset( -min( direction, zero3 ) );
		// odd slices start one slice "down"
		int oddSlicesStartIndex = evenSlicesStartIndex + gridOffsets.z;

		// set z to 0
		// adapt the local dimensions in such a way as to avoid out of bounds accesses at the "far corners"
		// the positive components of localSliceDirection affect the max corner of the slice
		// the negative ones the min corner (see *StartIndex). dimensions = max - min => use abs to count both correctly.
		localDimensions -= abs( localDirection );

		// always move 2 slices in local z direction, so we hit either odd or even slices in one kernel call
		gridOffsets.z *= 2;

		// there are floor( dimZ / 2 ) odd slices
		int numOddSlices = localDimensions.z / 2;
		int numEvenSlices = localDimensions.z - numOddSlices;

		int numCellsInSlice = localDimensions.x * localDimensions.y;

		// do all even slices
		Kernel_calculatePairLJForces<<<numEvenSlices * numCellsInSlice,1>>>(
				_positions.devicePtr(), _forces.devicePtr(),_cellInfos.devicePtr(), _domainValues.devicePtr(),
				evenSlicesStartIndex, make_int2( localDimensions ), gridOffsets,
				neighborOffset,
				epsilon, sigmaSquared, cutOffRadiusSquared
			);

		// do all odd slices
		Kernel_calculatePairLJForces<<<numOddSlices * numCellsInSlice,1>>>(
				_positions.devicePtr(), _forces.devicePtr(),_cellInfos.devicePtr(), _domainValues.devicePtr(),
				oddSlicesStartIndex, make_int2( localDimensions ), gridOffsets,
				neighborOffset,
				epsilon, sigmaSquared, cutOffRadiusSquared
			);
	}

	cellPairsTimer.end();

	singleCellsTimer.printElapsedTime( "intra cell LJ forces: %f ms " );
	cellPairsTimer.printElapsedTime( "inter cell LJ forces: %f ms\n" );
}

void LinkedCellsCUDA_Internal::reducePotentialAndVirial( OUT float &potential, OUT float &virial ) {
	const std::vector<unsigned long> &innerCellIndices = _linkedCells.getInnerCellIndices();
	const std::vector<unsigned long> &boundaryCellIndices = _linkedCells.getBoundaryCellIndices();

	potential = 0.0f;
	for( int i = 0 ; i < innerCellIndices.size() ; i++ ) {
		potential += _domainValues[ i ].x;
		virial += _domainValues[ i ].y;
	}
	for( int i = 0 ; i < boundaryCellIndices.size() ; i++ ) {
		potential += _domainValues[ i ].x;
		virial += _domainValues[ i ].y;
	}

	// every contribution is added twice so divide by 2
	potential /= 2.0f;
	virial /= 2.0f;

	// TODO: I have no idea why the sign is different in the GPU code...
	virial = -virial;
}
