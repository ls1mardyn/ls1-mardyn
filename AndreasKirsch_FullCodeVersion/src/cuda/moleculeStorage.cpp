/*
 * moleculeStorage.cpp
 *
 *  Created on: Jun 2, 2011
 *      Author: andreas
 */

#include "config.h"

#include "Domain.h"

#include "moleculeStorage.h"

#include "molecules/Molecule.h"

#include "math.h"

#include "benchmark.h"

void MoleculeStorage::uploadState() {
	_moleculePositions.clear();
	_moleculeQuaternions.clear();
#ifdef TEST_QUATERNION_MATRIX_CONVERSION
	_moleculeRotations.clear();
#endif

	_moleculeForces.clear();
	_moleculeTorque.clear();

	_moleculeComponentTypes.clear();
	_cellStartIndices.clear();

	int numCells = _linkedCells.getCells().size();

	unsigned numMolecules = 0;

	int maxNumMoleculesPerCell = 0;

	for( int i = 0 ; i < numCells ; i++ ) {
		const Cell &cell = _linkedCells.getCells()[i];

		_cellStartIndices.push( numMolecules );

		const std::list<Molecule*> &particles = cell.getParticlePointers();
#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
		for( int componentType = 0 ; componentType < _domain.getComponents().size() ; componentType++ ) {
#endif
			for( std::list<Molecule*>::const_iterator iterator = particles.begin() ; iterator != particles.end() ; iterator++ ) {
				Molecule &molecule = **iterator;

#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
				if( molecule.componentid() != componentType ) {
					continue;
				}
#endif

				_moleculeComponentTypes.push( molecule.componentid() );

				const floatType3 position = make_floatType3( molecule.r(0), molecule.r(1), molecule.r(2) );
				_moleculePositions.push( position );

				const Quaternion &dQuaternion = molecule.q();
				QuaternionStorage quaternion;
				quaternion.w = dQuaternion.qw();
				quaternion.x = dQuaternion.qx();
				quaternion.y = dQuaternion.qy();
				quaternion.z = dQuaternion.qz();
				_moleculeQuaternions.push( quaternion );

#ifdef TEST_QUATERNION_MATRIX_CONVERSION
				{
					Matrix3x3Storage rot;

					const floatType ww=quaternion.w*quaternion.w;
					const floatType xx=quaternion.x*quaternion.x;
					const floatType yy=quaternion.y*quaternion.y;
					const floatType zz=quaternion.z*quaternion.z;
					const floatType xy=quaternion.x*quaternion.y;
					const floatType zw=quaternion.z*quaternion.w;
					const floatType xz=quaternion.x*quaternion.z;
					const floatType yw=quaternion.y*quaternion.w;
					const floatType yz=quaternion.y*quaternion.z;
					const floatType xw=quaternion.x*quaternion.w;

					rot.rows[0].x=ww+xx-yy-zz;
					rot.rows[0].y=2*(xy-zw);
					rot.rows[0].z=2*(xz+yw);

					rot.rows[1].x=2*(xy+zw);
					rot.rows[1].y=ww-xx+yy-zz;
					rot.rows[1].z=2*(yz-xw);

					rot.rows[2].x=2*(xz-yw);
					rot.rows[2].y=2*(yz+xw);
					rot.rows[2].z=ww-xx-yy+zz;

					_moleculeRotations.push(rot);
				}
#endif

				numMolecules++;
			}
#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
		}
#endif

		maxNumMoleculesPerCell = std::max<int>( maxNumMoleculesPerCell, numMolecules - _cellStartIndices.getHostBuffer().back() );
	}

	printf( "average molecules per cell: %f (= %i / %i); max molecules per cell: %i\n", (float) numMolecules / numCells, numMolecules, numCells, maxNumMoleculesPerCell );

	_cellStartIndices.push( numMolecules );

	_moleculePositions.updateDevice();
	_moleculeQuaternions.updateDevice();

#ifndef TEST_QUATERNION_MATRIX_CONVERSION
	_moleculeRotations.resize( numMolecules );
#else
#	warning CPU testing quaternion matrix conversion
	_moleculeRotations.updateDevice();
#endif

	_moleculeForces.resize( numMolecules );
	_moleculeForces.zeroDevice();

	_moleculeTorque.resize( numMolecules );
	_moleculeTorque.zeroDevice();

	_moleculeComponentTypes.updateDevice();
	_cellStartIndices.updateDevice();

	_convertQuaternionsToRotations.call().
			parameter(numMolecules).
			setBlockShape(QUATERNION_BLOCK_SIZE, 1, 1).
			executeAtLeast((numMolecules + QUATERNION_BLOCK_SIZE - 1) / QUATERNION_BLOCK_SIZE);
}

struct CPUCudaVectorErrorMeasure {
	double totalCPUMagnitude, totalCudaMagnitude, totalSquareDeviation;
	double totalError, totalRelativeError;
	int numDataPoints;

	const char *name;

	CPUCudaVectorErrorMeasure(const char *name)
	: name(name), totalCPUMagnitude(0.0f), totalSquareDeviation(0.0f), totalCudaMagnitude(0.0f), totalError(0.0f), totalRelativeError(0.0f), numDataPoints(0) {
	}

	void registerErrorFor( const floatType3 &cpuResult, const floatType3 &cudaResult ) {
		const double epsilon = 5.96e-09f;

		// TODO: add convert_double3 macro/define!
#ifndef CUDA_DOUBLE_MODE
		const double3 delta = make_double3(cpuResult) - make_double3(cudaResult);
#else
		const double3 delta = cpuResult - cudaResult;
#endif

		const double cpuLength = length(cpuResult);
		const double cudaLength = length(cudaResult);
		const double deltaLengthSquared = dot(delta, delta);
		const double deltaLength = sqrt( deltaLengthSquared );

		totalCPUMagnitude += cpuLength;
		totalCudaMagnitude += cudaLength;
		totalSquareDeviation += deltaLengthSquared;

		totalError += deltaLength;
		if( cpuLength > epsilon ) {
			totalRelativeError += deltaLength / cpuLength;
		}
		numDataPoints++;
	}

	double getAverageCPUMagnitude() {
		return totalCPUMagnitude / numDataPoints;
	}

	double getAverageCUDAMagnitude() {
		return totalCudaMagnitude / numDataPoints;
	}

	double getAbsoluteAverageError() {
		return totalError / numDataPoints;
	}

	double getRelativeAverageError() {
		return totalRelativeError / numDataPoints;
	}

	double getRMSError() {
		return sqrt( numDataPoints * totalSquareDeviation ) / totalCPUMagnitude;
	}

	void report() {
		printf( "%s:\n"
				"  average CPU: %f\n"
				"  average CUDA: %f\n"
				"\n"
				"  average error: %f; average relative error: %f\n"
				"  RMS error: %f\n",
				name,
				getAverageCPUMagnitude(), getAverageCUDAMagnitude(),
				getAbsoluteAverageError(), getRelativeAverageError(),
				getRMSError()
			);
	}
};

void MoleculeStorage::compareResultsToCPURef( const std::vector<floatType3> &forces, const std::vector<floatType3> &torque ) {
	CPUCudaVectorErrorMeasure forceErrorMeasure( "force statistics" ), torqueErrorMeasure( "torque statistics" );

	const int numCells = _linkedCells.getCells().size();

	int currentIndex = 0;
	for( int i = 0 ; i < numCells ; i++ ) {
		const Cell &cell = _linkedCells.getCells()[i];

		const std::list<Molecule*> &particles = cell.getParticlePointers();

#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
		for( int componentType = 0 ; componentType < _domain.getComponents().size() ; componentType++ ) {
#endif
			for( std::list<Molecule*>::const_iterator iterator = particles.begin() ; iterator != particles.end() ; iterator++ ) {
				Molecule &molecule = **iterator;

#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
				if( molecule.componentid() != componentType ) {
					continue;
				}
#endif

				const floatType3 &cudaForce = forces[currentIndex];
				const floatType3 &cudaTorque = torque[currentIndex];
				currentIndex++;

				if( !cell.isHaloCell() ) {
					// we are going to compare F and M, so combine the sites
					molecule.calcFM();

					const floatType3 cpuForce = make_floatType3( molecule.F(0), molecule.F(1), molecule.F(2) );
					const floatType3 cpuTorque = make_floatType3( molecule.M(0), molecule.M(1), molecule.M(2) );

					forceErrorMeasure.registerErrorFor( cpuForce, cudaForce );
					torqueErrorMeasure.registerErrorFor( cpuTorque, cudaTorque );
				}

				// clear the molecule after comparing the values to make sure that only the GPU values are applied
				molecule.clearFM();
			}
#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
		}
#endif
	}

	forceErrorMeasure.report();
	torqueErrorMeasure.report();

	simulationStats.forceRMSError.addDataPoint( forceErrorMeasure.getRMSError() );
	simulationStats.torqueRMSError.addDataPoint( torqueErrorMeasure.getRMSError() );
}

void MoleculeStorage::downloadResults() {
#ifndef CUDA_UNPACKED_STORAGE
	const std::vector<floatType3> &forces = _moleculeForces.copyToHost();
	const std::vector<floatType3> &torque = _moleculeTorque.copyToHost();
#else
	std::vector<floatType3> forces;
	std::vector<floatType3> torque;

	_moleculeForces.copyToHost( forces );
	_moleculeTorque.copyToHost( torque );
#endif

#ifdef COMPARE_TO_CPU
	compareResultsToCPURef( forces, torque );
#endif

	const int numCells = _linkedCells.getCells().size();

	int currentIndex = 0;
	for( int i = 0 ; i < numCells ; i++ ) {
		const Cell &cell = _linkedCells.getCells()[i];

		const std::list<Molecule*> &particles = cell.getParticlePointers();
#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
		for( int componentType = 0 ; componentType < _domain.getComponents().size() ; componentType++ ) {
#endif
			for( std::list<Molecule*>::const_iterator iterator = particles.begin() ; iterator != particles.end() ; iterator++ ) {
				Molecule &molecule = **iterator;

#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
				if( molecule.componentid() != componentType ) {
					continue;
				}
#endif

				molecule.setF((floatType*) &forces[currentIndex]);
				molecule.setM((floatType*) &torque[currentIndex]);
				currentIndex++;
			}
#ifdef CUDA_SORT_CELLS_BY_COMPONENTTYPE
		}
#endif
	}
}
