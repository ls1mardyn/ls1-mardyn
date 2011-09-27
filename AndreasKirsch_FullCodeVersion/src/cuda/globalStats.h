/*
 * cellStats.h
 *
 *  Created on: Jun 6, 2011
 *      Author: andreas
 */

#ifndef GLOBALSTATS_H_
#define GLOBALSTATS_H_

#include "config.h"
#include "cudaComponent.h"

class GlobalStats : public CUDAInteractionCalculationComponent {
public:
	GlobalStats( const CUDAComponent &component )
		: CUDAInteractionCalculationComponent(component),
		  _cellStats( _module.getGlobal<CellStatsStorage *>("cellStats") ),
#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
		  _cellLocks( _module.getGlobal<LockStorage *>("cellLocks") ),
#endif
		  _potential( 0.0f ), _virial( 0.0f ) {
	}

	virtual void preInteractionCalculation();
	virtual void postInteractionCalculation();

	floatType getPotential() const {
		return _potential;
	}

	floatType getVirial() const {
		return _virial;
	}

protected:
	floatType _potential;
	floatType _virial;

	CUDA::PackedGlobalVector<CellStatsStorage> _cellStats;

#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
	CUDA::GlobalDeviceBuffer<LockStorage> _cellLocks;
#endif
};


#endif /* GLOBALSTATS_H_ */
