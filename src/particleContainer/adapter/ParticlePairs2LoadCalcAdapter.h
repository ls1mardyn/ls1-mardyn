#ifndef PARTICLEPAIRS2LOADCALCADAPTER_H_
#define PARTICLEPAIRS2LOADCALCADAPTER_H_

#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "utils/Logger.h"

//! @brief used for guessing load
//! @author Martin Buchholz
//!
class ParticlePairs2LoadCalcAdapter : public ParticlePairsHandler {
public:
	//! Constructor
	ParticlePairs2LoadCalcAdapter() {
	}

	ParticlePairs2LoadCalcAdapter(int globalCellsPerDim[3], int lowCorner[3], double cellSize[3], ParticleContainer* moleculeContainer) {
		for (int dim = 0; dim < 3; dim++) {
			_globalCellsPerDim[dim] = globalCellsPerDim[dim];
			_lowCorner[dim] = lowCorner[dim];
			_cellSize[dim] = cellSize[dim];
		}
		_globalNumCells = _globalCellsPerDim[0] * _globalCellsPerDim[1] * _globalCellsPerDim[2];
		_globalLoadPerCell = new float[_globalNumCells];

		for (int dim = 0; dim < 3; dim++) {
			_bBMin[dim] = moleculeContainer->getBoundingBoxMin(dim);
		}
		for (unsigned int i = 0; i < _globalNumCells; i++) {
			_globalLoadPerCell[i] = 0.0;
		}
		_localLoad = 0.0;
	}

	//! Destructor
	~ParticlePairs2LoadCalcAdapter() {
		delete[] _globalLoadPerCell;
	}

	//! @brief initialize
	void init() {
	}

	//! @brief finish
	void finish() {
		float* temp = new float[_globalNumCells];
		MPI_CHECK( MPI_Allreduce(_globalLoadPerCell, temp, _globalNumCells, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD) );
		delete[] _globalLoadPerCell;
		_globalLoadPerCell = temp;
		//cout << "LocalLoad: " << _localLoad << std::endl;
	}

	//! @brief calculate force between pairs and collect macroscopic contribution
	//!
	//! For all pairs, the force between the two Molecules has to be calculated
	//! and stored in the molecules. For original pairs(pairType 0), the contributions
	//! to the macroscopic values have to be collected
	double processPair(Molecule& particle1, Molecule& particle2, double distanceVector[3],
			PairType pairType, double dd, bool calculateLJ, double* force = NULL) {
		if (pairType == MOLECULE_MOLECULE) {
			int cellIndex[3]; // 3D Cell index (local)
			int globalCellIdx[3]; // 3D Cell index (global)

			for (int dim = 0; dim < 3; dim++) {
				cellIndex[dim] = (int) floor((particle1.r(dim) - _bBMin[dim]) / _cellSize[dim]);
				globalCellIdx[dim] = _lowCorner[dim] + cellIndex[dim];
			}
			unsigned long cellid = _globalCellsPerDim[0] * (globalCellIdx[2] * _globalCellsPerDim[1] + globalCellIdx[1]) + globalCellIdx[0];
			_globalLoadPerCell[cellid] += 1.0;
			_localLoad += 1.0;
		}
		return 0.0;
	}

	double processPair(Molecule& particle1, Molecule& particle2, double distanceVector[3], PairType pairType, double dd, bool calculateLJ)
	{
		return this->processPair(particle1, particle2, distanceVector, pairType, dd, calculateLJ, NULL);
	}

//	void recordRDF() {
//		return;
//	}

	void giveStatus() {
		std::cout << "Adapter: ParticlePairs2LoadCalcAdapter" << std::endl;
	}

	float* getLoad() {
		return _globalLoadPerCell;
	}

private:
	int _globalCellsPerDim[3];
	unsigned long _globalNumCells;
	float* _globalLoadPerCell;
	KDDecomposition* _kddecomp;
	double _bBMin[3];
	int _lowCorner[3];
	double _cellSize[3];
	float _localLoad;
};

#endif /*PARTICLEPAIRS2LOADCALCADAPTER_H_*/
