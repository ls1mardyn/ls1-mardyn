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
		  _cellLocks( _module.getGlobal<LockStorage *>("cellLocks") ),
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

	CUDA::GlobalDeviceBuffer<LockStorage> _cellLocks;
};


#endif /* GLOBALSTATS_H_ */
