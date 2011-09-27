/*
 * moleculeStorage.h
 *
 *  Created on: May 28, 2011
 *      Author: andreas
 */

#ifndef MOLECULESTORAGE_H_
#define MOLECULESTORAGE_H_

#include "cudaComponent.h"

#include "config.h"

#include "sharedDecls.h"

class MoleculeStorage : public CUDAInteractionCalculationComponent {
public:
	MoleculeStorage( const CUDAComponent &component ) :
		CUDAInteractionCalculationComponent( component ),
#ifndef CUDA_UNPACKED_STORAGE
		_moleculePositions( _module.getGlobal<floatType3 *>("moleculePositions") ),
		_moleculeQuaternions( _module.getGlobal<QuaternionStorage *>("moleculeQuaternions") ),
		_moleculeRotations( _module.getGlobal<Matrix3x3Storage *>("moleculeRotations") ),

		_moleculeForces( _module.getGlobal<floatType3 *>("moleculeForces") ),
		_moleculeTorque( _module.getGlobal<floatType3 *>("moleculeTorque") ),
#else
		_moleculePositions( _module.getGlobalArray<floatType *, 3>("moleculePositions") ),
		_moleculeQuaternions( _module.getGlobalArray<floatType *, 4>("moleculeQuaternions") ),
		_moleculeRotations( _module.getGlobalArray<floatType *, 9>("moleculeRotations") ),

		_moleculeForces( _module.getGlobalArray<floatType *, 3>("moleculeForces") ),
		_moleculeTorque( _module.getGlobalArray<floatType *, 3>("moleculeTorque") ),
#endif

		_moleculeComponentTypes( _module.getGlobal<Molecule_ComponentType *>("moleculeComponentTypes") ),
		_cellStartIndices( _module.getGlobal<unsigned *>("cellStartIndices") ),

		_convertQuaternionsToRotations( _module.getFunction( "convertQuaternionsToRotations" ) )
		{
	}

	void preInteractionCalculation() {
		uploadState();
	}

	void postInteractionCalculation() {
		downloadResults();
	}

protected:
	void uploadState();
	void downloadResults();

	void compareResultsToCPURef( const std::vector<floatType3> &forces, const std::vector<floatType3> &torque );

#ifndef CUDA_UNPACKED_STORAGE
	CUDA::PackedGlobalVector<floatType3> _moleculePositions, _moleculeForces, _moleculeTorque;
	CUDA::PackedGlobalVector<QuaternionStorage> _moleculeQuaternions;
	CUDA::PackedGlobalVector<Matrix3x3Storage> _moleculeRotations;
#else
	CUDA::UnpackedGlobalVector<floatType, 3> _moleculePositions, _moleculeForces, _moleculeTorque;
	CUDA::UnpackedGlobalVector<floatType, 4> _moleculeQuaternions;

	CUDA::UnpackedGlobalVector<floatType, 9> _moleculeRotations;
#endif

	CUDA::PackedGlobalVector<Molecule_ComponentType> _moleculeComponentTypes;
	CUDA::PackedGlobalVector<unsigned> _cellStartIndices;

	CUDA::Function _convertQuaternionsToRotations;
};

#endif /* MOLECULESTORAGE_H_ */
