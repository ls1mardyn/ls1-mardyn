/*
 * moleculeInteraction.h
 *
 *  Created on: Jun 8, 2011
 *      Author: andreas
 */

#ifndef MOLECULEINTERACTION_H_
#define MOLECULEINTERACTION_H_

#include "helpers.h"

#include "cudaComponent.h"
#include "componentDescriptor.h"
#include "globalStats.h"
#include "moleculeStorage.h"
#include "moleculePairHandler.h"
#include "domainTraverser.h"
#include "config.h"

class MoleculeInteraction : public CUDAComponent {
public:
	MoleculeInteraction( const CUDA::Module &module, LinkedCells &linkedCells, Domain &domain ) :
		CUDAComponent(module, linkedCells, domain),

		_globalStats( *this ),
		_moleculeStorage( *this ),

		_componentDescriptorStorage( *this ),
		_moleculePairHandler( *this ),

		_domainTraverser( *this ),

#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
		_createSchedulers( module.getFunction("createSchedulers") ),
		_destroySchedulers( module.getFunction("destroySchedulers") ),
#endif

		_cellPairProcessor( module.getFunction("processCellPair") ),
		_cellProcessor( module.getFunction("processCell") )

	{
		// upload data from CUDAStaticDataComponents
		_moleculePairHandler.upload();
		_componentDescriptorStorage.upload();
	}

	void calculate(double &potential, double &virial) {
		CUDAFrameTimer.begin();

		// pre-force calculation handling from CUDAForceCalculationComponents
		CUDAPreTimer.begin();
		_globalStats.preInteractionCalculation();
		_moleculeStorage.preInteractionCalculation();
		CUDAPreTimer.end();

		CUDATotalProcessingTimer.begin();
		CUDAPairProcessingTimer.begin();
		for( int stageIndex = 0 ; stageIndex < _domainTraverser.getInterCellStageCount() ; stageIndex++ ) {
			_domainTraverser.preInterCellStage( stageIndex );

			for( int subStageIndex = 0 ; subStageIndex < _domainTraverser.getInterCellSubStageCount( stageIndex ) ; subStageIndex++ ) {
				printf( "inter-cell: %i:%i\n", stageIndex, subStageIndex );
				_domainTraverser.preInterCellSubStage( stageIndex, subStageIndex );

#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
				_createSchedulers.call().execute();
#endif
				_cellPairProcessor.call().
					setBlockShape( WARP_SIZE, NUM_WARPS, 1 ).
#ifndef CUDA_WARP_BLOCK_CELL_PROCESSOR
					executeAtLeast( _domainTraverser.getInterCellJobCount(stageIndex, subStageIndex) );
#else
					execute( 16 );
#endif

#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
				_destroySchedulers.call().execute();
#endif
			}
		}
		CUDAPairProcessingTimer.end();

#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
		_createSchedulers.call().execute();
#endif

		CUDASingleProcessingTimer.begin();
		printf( "intra-cell\n" );

		_domainTraverser.preIntraCellStage();
		_cellProcessor.call().
				setBlockShape( WARP_SIZE, NUM_WARPS, 1 ).
#ifndef CUDA_WARP_BLOCK_CELL_PROCESSOR
				executeAtLeast( _domainTraverser.getIntraCellJobCount() );
#else
				execute( 16 );
#endif

		CUDASingleProcessingTimer.end();

#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
		_destroySchedulers.call().execute();
#endif

		CUDATotalProcessingTimer.end();

		// post-force calculation handling from CUDAForceCalculationComponents
		CUDAPostTimer.begin();
		_globalStats.postInteractionCalculation();
		_moleculeStorage.postInteractionCalculation();
		CUDAPostTimer.end();

		potential = _globalStats.getPotential();
		virial = _globalStats.getVirial();

		CUDAFrameTimer.end();

		simulationStats.CUDA_frameTime.addDataPoint(CUDAFrameTimer.getElapsedTime());
		simulationStats.CUDA_preTime.addDataPoint(CUDAPreTimer.getElapsedTime());
		simulationStats.CUDA_postTime.addDataPoint(CUDAPostTimer.getElapsedTime());
		simulationStats.CUDA_pairTime.addDataPoint(CUDAPairProcessingTimer.getElapsedTime());
		simulationStats.CUDA_singleTime.addDataPoint(CUDASingleProcessingTimer.getElapsedTime());
		simulationStats.CUDA_processingTime.addDataPoint(CUDATotalProcessingTimer.getElapsedTime());
	}

protected:
	CUDA::Function _cellPairProcessor, _cellProcessor;
#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
	CUDA::Function _createSchedulers, _destroySchedulers;
#endif

	GlobalStats _globalStats;
	MoleculeStorage _moleculeStorage;
	MoleculePairHandler _moleculePairHandler;
	ComponentDescriptorStorage _componentDescriptorStorage;
	DomainTraverser _domainTraverser;

	CUDA::EventTimer CUDAFrameTimer, CUDAPreTimer, CUDATotalProcessingTimer, CUDAPairProcessingTimer, CUDASingleProcessingTimer, CUDAPostTimer;
};

#endif /* MOLECULEINTERACTION_H_ */
