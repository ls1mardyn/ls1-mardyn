#ifndef CUDA_PAIRTRAVERSER_H_
#define CUDA_PAIRTRAVERSER_H_

// we want to use the CUDA types: int3 etc
#include <vector_types.h>
#include "cutil_math.h"

#include "cudaComponent.h"

class DomainTraverser : public CUDAComponent {
public:
	DomainTraverser( const CUDAComponent &component )
		:
			CUDAComponent(component),

			_startIndex( _module.getGlobal<int>("_ZN15DomainTraverser10startIndexE") ),
			_dimension( _module.getGlobal<int2>("_ZN15DomainTraverser9dimensionE") ),
			_gridOffsets( _module.getGlobal<int3>("_ZN15DomainTraverser11gridOffsetsE") ),
			_neighborOffset( _module.getGlobal<int>("_ZN15DomainTraverser14neighborOffsetE") ),

			_numJobs( _module.getGlobal<int>("_ZN15DomainTraverser7numJobsE") )
	{
		initStages();
	}

	int getInterCellStageCount() const {
		return numDirections;
	}

	int getInterCellSubStageCount(int stageIndex) const {
		return 2;
	}

	void preInterCellStage(int stageIndex) {
		const InterCellStage &interCellStage = interCellStages[stageIndex];

		_gridOffsets.set( interCellStage.gridOffsets );
		_dimension.set( interCellStage.localDimensions );
		_neighborOffset.set( interCellStage.neighborOffset );
	}

	void preInterCellSubStage(int stageIndex, int subStageIndex) {
		const InterCellStage::SubStage &subStage = interCellStages[stageIndex].subStages[subStageIndex];
		_numJobs.set( subStage.numCellPairs );
		_startIndex.set( subStage.startIndex );
	}

	void preIntraCellStage() {
		_gridOffsets.set( intraCellStage.gridOffsets );
		_dimension.set( intraCellStage.localDimensions );

		_numJobs.set( intraCellStage.numCells );
		_startIndex.set( intraCellStage.startIndex );
	}

	int getInterCellJobCount(int stageIndex, int subStageIndex) const {
		return interCellStages[stageIndex].subStages[subStageIndex].numCellPairs;
	}

	int getIntraCellJobCount() {
		return intraCellStage.numCells;
	}

protected:
	static const int numDirections = 13;

	struct InterCellStage {
		struct SubStage {
			int numCellPairs;
			int startIndex;
		} subStages[2];

		int2 localDimensions;
		int3 gridOffsets;
		int neighborOffset;
	};

	InterCellStage interCellStages[numDirections];

	struct IntraCellStage {
		int numCells;
		int startIndex;

		int2 localDimensions;
		int3 gridOffsets;
	};

	IntraCellStage intraCellStage;

	int getDirectionOffset( const int3 &direction ) {
		return _linkedCells.cellIndexOf3DIndex( direction.x, direction.y, direction.z );
	}

	int getCellOffset( const int3 &cell ) {
		return getDirectionOffset( cell );
	}

	void initStages() {
		const int3 dimensions = *(int3*) _linkedCells.getCellDimensions();
		assert( dimensions.x > 2 && dimensions.y > 2 && dimensions.z > 2 );

		const int3 haloWidth = make_int3( 1, 1, 1 );

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

		const int3 cellDirections = {
				getDirectionOffset( xDirection ),
				getDirectionOffset( yDirection ),
				getDirectionOffset( zDirection )
			};

		intraCellStage.gridOffsets = cellDirections;
		intraCellStage.localDimensions = make_int2( dimensions );

		intraCellStage.startIndex = getCellOffset( haloWidth );
		intraCellStage.numCells = (dimensions.x - 2) * (dimensions.y - 2) * (dimensions.z - 2);

		for( int i = 0 ; i < numDirections ; i++ ) {
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
				localDimensions = make_int3( dimensions.y, dimensions.z, dimensions.x );
				gridOffsets = make_int3( cellDirections.y, cellDirections.z, cellDirections.x );
			}
			else if( direction.y == 1 ) {
				// xz plane (main direction: y)
				localDirection = make_int3( direction.x, direction.z, direction.y );
				localDimensions = make_int3( dimensions.x, dimensions.z, dimensions.y );
				gridOffsets = make_int3( cellDirections.x, cellDirections.z, cellDirections.y );
			}
			else if( direction.z == 1 ) {
				// xy plane (main direction: z)
				localDirection = direction;
				localDimensions = dimensions;
				gridOffsets = cellDirections;
			}
			else {
				assert( false );
			}

			// determine the startOffset as first cell near (0,0,0)
			// so that start + neighborOffset won't be out of bounds
			int evenSlicesStartIndex = getCellOffset( -min( direction, zero3 ) );

			// odd slices start one slice "down"
			int oddSlicesStartIndex = evenSlicesStartIndex + gridOffsets.z;

			// adapt the local dimensions in such a way as to avoid out of bounds accesses at the "far corners"
			// the positive components of localSliceDirection affect the max corner of the slice
			// the negative ones the min corner (see *StartIndex). dimensions = max - min => use abs to count both correctly.
			localDimensions -= abs( localDirection );

			// always move 2 slices in local z direction, so we hit either odd or even slices in one kernel call
			gridOffsets.z *= 2;

			// there are floor( dimZ / 2 ) odd slices
			const int numOddSlices = localDimensions.z / 2;
			const int numEvenSlices = localDimensions.z - numOddSlices;

			const int numCellsInSlice = localDimensions.x * localDimensions.y;

			// set the stage parameters
			InterCellStage &interCellStage = interCellStages[i];

			interCellStage.gridOffsets = gridOffsets;
			interCellStage.localDimensions = make_int2( localDimensions );
			interCellStage.neighborOffset = neighborOffset;

			// even slices
			InterCellStage::SubStage &evenSubStage = interCellStage.subStages[0];
			evenSubStage.numCellPairs = numEvenSlices * numCellsInSlice;
			evenSubStage.startIndex = evenSlicesStartIndex;

			// odd slices
			InterCellStage::SubStage &oddSubStage = interCellStage.subStages[1];
			oddSubStage.numCellPairs = numOddSlices * numCellsInSlice;
			oddSubStage.startIndex = oddSlicesStartIndex;
		}
	}

	void testStages() const {
		struct InteractionLog {
			int numInteractions;
			int directionMask;

			InteractionLog() : numInteractions( 0 ), directionMask( 0 ) {}
		};
		int numCells = _linkedCells.getCells().size();
		InteractionLog *interactionLog = new InteractionLog[ numCells ];

		for( int stage = 0 ; stage < getInterCellStageCount() ; stage++ ) {
			for( int subStage = 0 ; subStage < getInterCellSubStageCount( stage ) ; subStage++ ) {
				for( int job = 0 ; job < getInterCellJobCount( stage, subStage ) ; job++ ) {
					const InterCellStage &params = interCellStages[stage];
					const int3 gridIndex = make_int3( job % params.localDimensions.x,
								(job / params.localDimensions.x) % params.localDimensions.y,
								job / params.localDimensions.x / params.localDimensions.y
							);

					const int cellIndex = params.subStages[subStage].startIndex + dot( gridIndex, params.gridOffsets );
					const int neighborIndex = cellIndex + params.neighborOffset;

					if( cellIndex < 0 || cellIndex >= numCells ) {
						printf( "Bad cell index: %i in stage/subStage/job: %i/%i/%i\n", cellIndex, stage, subStage, job );
						continue;
					}
					if( neighborIndex < 0 || neighborIndex >= numCells ) {
						printf( "Bad neighbor index: %i in stage/subStage/job: %i/%i/%i\n", neighborIndex, stage, subStage, job );
						continue;
					}

					InteractionLog &A = interactionLog[ cellIndex ];
					InteractionLog &B = interactionLog[ neighborIndex ];

					A.numInteractions++;
					B.numInteractions++;
					const int directionBit = 1 << (stage * 2 + subStage);
					A.directionMask |= directionBit;
					B.directionMask |= directionBit;
				}
			}
		}

		const std::vector<unsigned long> &innerCellIndices = _linkedCells.getInnerCellIndices();
		const std::vector<unsigned long> &boundaryCellIndices = _linkedCells.getBoundaryCellIndices();

		for( int i = 0 ; i < innerCellIndices.size() ; i++ ) {
			const int cellIndex = innerCellIndices[i];
			const InteractionLog &cell = interactionLog[ cellIndex ];
			if( cell.numInteractions != 26 ) {
				int x, y, z;
				_linkedCells.cellIndexTo3DIndex( cellIndex, x, y, z );
				printf( "interaction test failed for cell %i (%i, %i, %i): numInteractions = %i, directionMask = 0x%x \n", cellIndex, x, y, z, cell.numInteractions, cell.directionMask );
			}
		}
		for( int i = 0 ; i < boundaryCellIndices.size() ; i++ ) {
			const int cellIndex = boundaryCellIndices[i];
			const InteractionLog &cell = interactionLog[ cellIndex ];
			if( cell.numInteractions != 26 ) {
				int x, y, z;
				_linkedCells.cellIndexTo3DIndex( cellIndex, x, y, z );
				printf( "interaction test failed for cell %i (%i, %i, %i): numInteractions = %i, directionMask = 0x%x \n", cellIndex, x, y, z, cell.numInteractions, cell.directionMask );
			}
		}

		delete[] interactionLog;
	}

	CUDA::Global<int> _startIndex;
	CUDA::Global<int2> _dimension;
	CUDA::Global<int3> _gridOffsets;
	CUDA::Global<int> _neighborOffset;

	CUDA::Global<int> _numJobs;

};

#endif
